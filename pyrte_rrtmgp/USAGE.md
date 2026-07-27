# The Simple Spectral Model (SSM)

See [Williams (2026), *Bridging clarity and accuracy: A simple spectral
longwave radiation scheme for idealized climate modeling*](https://doi.org/10.1029/2025MS005405)
and Czarnecki and Pincus (2026) for the model formulation, and the Fortran
implementation in
[rte-rrtmgp `mo_optics_ssm.F90`](https://github.com/earth-system-radiation/rte-rrtmgp)
for the reference version this Python port reproduces.

## The model in two equations

Each absorbing feature (a *tag*) contributes a mass absorption coefficient
that decays exponentially away from a central wavenumber, so that
log&nbsp;&kappa; is triangle-shaped in wavenumber space:

```
kappa(nu) = kappa0 * exp(-|nu - nu0| / l)
```

The optical depth of each layer is the pressure-scaled sum over tags:

```
tau(nu) = (p_layer / pref) * sum_tags[ layer_mass * kappa(nu) ]
```

Planck sources are `B_nu(T) * dnus` evaluated on the same wavenumber grid, and
broadband fluxes are quadrature sums over that grid.

## Workflow

### 1. Configure the spectroscopy (`triangles`)

The spectroscopy is an `xarray.Dataset` holding one variable, `triangles`,
with dimensions `("tags", "params")`. The three params configure each feature:

| param    | meaning                                     | units    | constraint |
|----------|---------------------------------------------|----------|------------|
| `nu0`    | central wavenumber of the feature           | cm⁻¹     | —          |
| `l`      | spectral decay length away from `nu0`       | cm⁻¹     | > 0        |
| `kappa0` | peak mass absorption coefficient at `pref`  | m² kg⁻¹  | ≥ 0        |

Tag names identify the gas: the text before the first hyphen is the species
(`"h2o-rot"` → `h2o`), which must be one of `h2o`, `co2`, or `o3`. A gas may
own several tags (for example separate rotational and vibrational bands).

Two ready-made spectroscopies ship with the module:

* {data}`~pyrte_rrtmgp.ssm.SSM_W26` — Williams (2026); identical to the
  longwave default (`triangle_params_def_lw`) of the Fortran RTE-SSM:

  | tag       | nu0 (cm⁻¹) | l (cm⁻¹) | kappa0 (m² kg⁻¹) |
  |-----------|-----------|----------|------------------|
  | `co2`     | 667       | 12       | 110              |
  | `h2o-rot` | 0         | 64       | 282              |
  | `h2o-vr`  | 1600      | 52       | 24               |

* {data}`~pyrte_rrtmgp.ssm.SSM_CP26` — Czarnecki and Pincus (2026), which
  adds an `h2o-cont` continuum tag.

To define your own spectroscopy, build a Dataset in the same layout:

```python
import numpy as np
import xarray as xr

my_ssm = xr.Dataset(
    coords={
        "tags": ["co2", "h2o-rot"],
        "params": ["nu0", "l", "kappa0"],
    },
    data_vars={
        "triangles": (
            ["tags", "params"],
            np.array(
                [
                    [667.0, 12.0, 110.0],  # co2
                    [0.0, 64.0, 282.0],  # h2o-rot
                ]
            ),
        )
    },
)
```

### 2. Choose a spectral grid (`nus`, `dnus`)

`nus` are the wavenumbers (cm⁻¹, strictly increasing) at which absorption and
the Planck function are evaluated; `dnus` is the spectral width (cm⁻¹,
positive) credited to each point. Together they control spectral resolution
and total band coverage. The Fortran RTE-SSM default uses 41 points from 50 to
3000 cm⁻¹, with band edges at the midpoints between points and the outer edges
at 0 and 3500 cm⁻¹:

```python
nus_v = np.linspace(50.0, 3000.0, 41)
mids = 0.5 * (nus_v[:-1] + nus_v[1:])
dnus_v = np.diff(np.concatenate([[0.0], mids, [3500.0]]))

nus = xr.DataArray(nus_v, dims="gpt")
dnus = xr.DataArray(dnus_v, dims="gpt")
```

### 3. Initialize {class}`~pyrte_rrtmgp.ssm.GasOptics`

```python
from pyrte_rrtmgp.ssm import GasOptics, SSM_W26

ssm = GasOptics(
    spectral_data=SSM_W26,
    nus=nus,
    dnus=dnus,
    pref=SSM_W26.attrs["pref"] * 100.0,  # hPa -> Pa
)
```

```{warning}
The constructor's `pref` argument is in **Pa**, but the `pref` attribute
stored on `SSM_W26` (500.0) and `SSM_CP26` (1000.22) is in **hPa**. Multiply
by 100 when passing it on, as above. Passing the attribute directly runs
without error but scales all optical depths by 1/100.
```

### 4. Define the atmosphere

The input Dataset uses the same variable names as
{class}`pyrte_rrtmgp.rrtmgp.GasOptics`:

* `pres_layer`, `pres_level` — pressure at layer centers and boundaries (Pa)
* `temp_layer`, `temp_level` — temperature at layer centers and boundaries (K)
* `surface_temperature` — skin temperature (K), no vertical dimension
* one volume mixing ratio variable per species used by the spectroscopy
  (e.g. `h2o`, `co2`), in mol/mol. Well-mixed gases may be a scalar per
  column; they are broadcast across layers.

Levels are layer boundaries, so the `level` dimension is one element longer
than `layer`. Pressure may increase or decrease with index; the orientation is
detected automatically and recorded in the `top_at_1` attribute of the output.

```{warning}
The vertical (`layer`/`level`) axis must be the **last** dimension of every
variable. For multi-column data put the column dimension first, e.g.
`atm = atm.transpose("col", ...)`.
```

### 5. Compute optics with `compute()`

```python
optics = ssm.compute(atm, add_to_input=False)
```

The result contains `tau`, `layer_source`, `level_source`, `surface_source`,
`surface_source_jacobian`, and the spectral grid (`nus`, `dnus`) — everything
the longwave solver needs. With `add_to_input=True` the fields are written
into the input Dataset in place and `None` is returned.

### 6. Solve the radiative transfer equations

Merge in a broadband `surface_emissivity` and call the usual RTE accessor:

```python
problem = xr.merge([optics, atm.surface_emissivity])
fluxes = problem.rte.solve(add_to_input=False)
```

`fluxes` contains `lw_flux_up` and `lw_flux_down` on levels, in W/m².

## Worked example

A complete, self-contained script — a two-column, 20-layer idealized
atmosphere with H₂O and CO₂, solved end to end:

```python
import numpy as np
import xarray as xr

import pyrte_rrtmgp.rte  # noqa: F401  (registers the .rte accessor)
from pyrte_rrtmgp.ssm import SSM_W26, GasOptics

# --- 1. Spectral grid: 41 points from 50 to 3000 cm-1, band edges at the
#        midpoints between points, with the outer edges at 0 and 3500 cm-1.
nus_v = np.linspace(50.0, 3000.0, 41)
mids = 0.5 * (nus_v[:-1] + nus_v[1:])
dnus_v = np.diff(np.concatenate([[0.0], mids, [3500.0]]))

nus = xr.DataArray(nus_v, dims="gpt")
dnus = xr.DataArray(dnus_v, dims="gpt")

# --- 2. Gas optics from the bundled Williams (2026) triangles.
#        SSM_W26.attrs["pref"] is in hPa; the constructor expects Pa.
ssm = GasOptics(
    spectral_data=SSM_W26,
    nus=nus,
    dnus=dnus,
    pref=SSM_W26.attrs["pref"] * 100.0,
)

# --- 3. A small synthetic atmosphere: 2 columns x 20 layers.
#        Column dimension first; the vertical axis must be last.
ncol, nlay = 2, 20
plev = np.linspace(1000e2, 10e2, nlay + 1)  # Pa, surface first
play = 0.5 * (plev[:-1] + plev[1:])

tsfc = np.array([300.0, 285.0])
tlev = tsfc[:, None] - 65.0 * (1.0 - plev / plev[0])[None, :]
tlay = 0.5 * (tlev[:, :-1] + tlev[:, 1:])

h2o = 0.02 * (play / play[0]) ** 3  # moist near the surface, dry aloft

atm = xr.Dataset(
    data_vars={
        "pres_level": (("col", "level"), np.tile(plev, (ncol, 1))),
        "pres_layer": (("col", "layer"), np.tile(play, (ncol, 1))),
        "temp_level": (("col", "level"), tlev),
        "temp_layer": (("col", "layer"), tlay),
        "h2o": (("col", "layer"), h2o[None, :] * np.array([[1.0], [0.5]])),
        "co2": (("col",), np.array([420e-6, 420e-6])),  # well-mixed scalar
        "surface_temperature": (("col",), tsfc),
        "surface_emissivity": (("col",), np.array([0.98, 0.98])),
    }
)

# --- 4. Optical depth + Planck sources, then solve the RTE.
optics = ssm.compute(atm, add_to_input=False)
problem = xr.merge([optics, atm.surface_emissivity])
fluxes = problem.rte.solve(add_to_input=False)

# --- 5. Pressure decreases with index here, so level 0 is the surface
#        and the last level is the top of atmosphere.
olr = fluxes.lw_flux_up.isel(level=-1)
sfc_dn = fluxes.lw_flux_down.isel(level=0)

for c in range(ncol):
    print(
        f"column {c}: OLR = {float(olr[c]):6.1f} W/m2, "
        f"surface downwelling LW = {float(sfc_dn[c]):6.1f} W/m2"
    )
```

Output:

```
column 0: OLR =  338.0 W/m2, surface downwelling LW =  312.9 W/m2
column 1: OLR =  278.2 W/m2, surface downwelling LW =  247.7 W/m2
```

The warm, moist column emits more to space and returns more longwave to the
surface than the cool column with half the water vapor — the expected
behavior.

## Validation notes

* This Python implementation reproduces the Fortran RTE-SSM reference fluxes
  to better than 0.04 W/m² (relative differences below 0.02%) across the
  RFMIP, RCE, and CKDMIP example atmospheres.
* Against RRTMGP (`LW_G256`), the SSM with the `SSM_W26` defaults shows net
  longwave biases of roughly 20–30 W/m² at the top of atmosphere and
  40–50 W/m² at the surface on the same atmospheres. The SSM trades this
  accuracy for transparency and speed; it is an idealized model, not a
  replacement for RRTMGP in production radiation budgets.

## Limitations

* Longwave only; absorption only (no scattering, no clouds, no aerosols).
* Supported absorbers: `h2o`, `co2`, `o3`.
* Spectroscopy validation requires unique tags, finite triangle parameters,
  `l > 0`, `kappa0 >= 0`, strictly increasing `nus` (at least two points),
  and positive `dnus` of matching length.
