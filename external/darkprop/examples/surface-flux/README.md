# Surface-Flux DarkProp Runs

This directory contains the DarkProp setup used to compare the Boltzmann-equation
iteration with Monte Carlo propagation for 5 GeV halo dark matter.

## Physics and Geometry Choices

- DM model: `dm_model = "SI"`.
- Input type: `type = "flux"`.
- Medium: modified `HomoEarth`.
- Detector depth: `depth = [2.4]` km.
- DM mass: `mchi.txt` contains `5 GeV`.
- Surface flux: `surface_flux.txt` is the benchmark halo surface spectrum converted to
  `dPhi/(dT dOmega)` in `GeV^-1 cm^-2 s^-1 sr^-1`.
- The surface flux is independent of `sigma_chiN`; for a new cross section, only
  the sigma file and DarkProp config need to change.

Local DarkProp source edits already made:

- `darkprop/HomoEarth.hpp`: removed the `/ 0.985` normalization from crust target
  densities to match the benchmark homogeneous-earth composition.
- `darkprop/HomoEarth.hpp`: changed Fe from `57*constants::mu` to
  `56*constants::mu`.

## Reference MC Run Outputs

The full Monte Carlo HDF5 files and exported tables are generated under `out/`
when the DarkProp runs are executed.

Reference runs used for comparison:

- `sigma_chiN = 1e-33 cm^2`
  - Total MC sample: `50,000,000` detector-sphere crossings.
  - Merged table:
    `out/surface-flux-si-compare/darkprop_mc_dphidv.dat`
  - HDF5 inputs:
    `out/surface-flux-si-compare/*.hdf5`,
    `out/surface-flux-si-compare-seed101/*.hdf5`,
    `out/surface-flux-si-compare-seed201-40m/*.hdf5`
  - Ignore `out/surface-flux-si-compare-seed2`: seed 2 overlaps with seed 1 under
    DarkProp's `random_seed + MPI rank` rule.

- `sigma_chiN = 1e-32 cm^2`
  - Total MC sample: `10,000,000` detector-sphere crossings.
  - Merged table:
    `out/surface-flux-si-sigma1e-32-10m/darkprop_mc_dphidv.dat`
  - HDF5 input:
    `out/surface-flux-si-sigma1e-32-10m/surface-flux-si-sigma1e-32-10m_mchi5.000e+00GeV_sigma1.000e-32cm2.hdf5`

Low-stat exploratory run:

- `sigma_chiN = 1e-31 cm^2`
  - Config: `config_surface_flux_si_sigma1e-31_50m.toml`
  - This hit `maximal_simulation = 100000000` before collecting many crossings:
    only `Nsample = 3804`.
  - The exported table has only a few nonzero bins and should not be used as a
    final comparison.

## Run a New Cross Section

Use PowerShell from the repository root:

```powershell
cd <repo-root>
```

Create a sigma file, for example `3e-32 cm^2`:

```powershell
Set-Content -Path "external\darkprop\examples\surface-flux\sigma_3e-32.txt" -Value "3.0000000000000000e-32"
```

Copy one of the existing configs and edit these fields:

- `id`: unique output directory name.
- `sigma_file`: point to the new sigma file.
- `sample_size`: target number of detector-sphere crossings.
- `maximal_simulation`: cap on injected particles. For stronger cross sections,
  this must be much larger than `sample_size`.
- `random_seed`: keep it separated by more than the MPI rank count from previous
  runs. With `mpiexec -n 8`, use seeds like `601`, `701`, `801`, not consecutive
  seeds.

Example config fields:

```toml
id = "surface-flux-si-sigma3e-32-10m"

[input]
dm_model = "SI"
medium_model = "HomoEarth"
type = "flux"
mchi_file = "./mchi.txt"
sigma_file = "./sigma_3e-32.txt"
Tcut = 1.7073777777777778e-13

[input.flux]
flux_file = "./surface_flux.txt"
point_source = false
Tmin = 1.7073777777777778e-13
Tmax = 1.7073777777777781e-05

[output]
outdir = "./out"
overwrite = false
maximal_simulation = 200000000
sample_size = 10000000
depth = [2.4]
chunk_size = 100000
full_tracks = 0

[output.log]
screen_log_lag = 1000000

[numerical]
random_seed = 601
```

Run DarkProp in WSL:

```powershell
$surfaceFluxDir = wsl wslpath -a (Resolve-Path "external\darkprop\examples\surface-flux")
wsl bash -lc "cd '$surfaceFluxDir' && mpiexec -n 8 ../../build/install/bin/darkprop config_surface_flux_si_sigma3e-32-10m.toml"
```

## Export `dPhi/dv`

Use the exporter from the repo root. It outputs columns:

```text
v_over_c dPhi_dv_cm^-2_s^-1 stat_error_cm^-2_s^-1
```

Single HDF5:

```powershell
python codes\export_darkprop_mc_dphidv.py `
  --darkprop external\darkprop\examples\surface-flux\out\surface-flux-si-sigma3e-32-10m\surface-flux-si-sigma3e-32-10m_mchi5.000e+00GeV_sigma3.000e-32cm2.hdf5 `
  --nbins 180 `
  --output external\darkprop\examples\surface-flux\out\surface-flux-si-sigma3e-32-10m\darkprop_mc_dphidv.dat
```

Multiple independent runs can be accumulated by listing all HDF5 files:

```powershell
python codes\export_darkprop_mc_dphidv.py `
  --darkprop run1.hdf5 run2.hdf5 run3.hdf5 `
  --nbins 180 `
  --output merged_darkprop_mc_dphidv.dat
```

The exporter uses `--min-abs-cos 1e-3` by default. This drops extremely tangent
detector-sphere crossings where the official surface reconstruction weight
contains `1/|cos(theta)|`. In the `1e-33` run it removed only 32 events out of
50M but removed large artificial spikes. Use `--min-abs-cos 0` to reproduce the
raw estimator.

## Regenerating the Surface Flux

Usually this is not needed. To regenerate the benchmark surface-flux input:

```powershell
python codes\generate_darkprop_surface_flux.py `
  --mchi-mev 5000 `
  --sigma-chin-cm2 1e-33 `
  --n-energy 400 `
  --v-min-fraction 1e-4
```

This rewrites `mchi.txt`, `sigma.txt`, and `surface_flux.txt`.
