# Boltzmann Equation Earth Attenuation Codes

This repository contains the numerical codes used for the Boltzmann-equation
calculation of dark-matter attenuation through the Earth.

The repository is a cleaned code release extracted from the working notes
directory. Generated arrays, plots, logs, LaTeX files, and local build products
are intentionally excluded.

## Repository layout

```text
codes/
  iteration_dirac_dm_isoscalar_vector_coupling.py
  run_dirac_dm_isoscalar_vector_coupling.py
  verify_iteration_integrals.py
  plot_detector_cumulative_spectra.py
  generate_darkprop_surface_flux.py
  export_darkprop_mc_dphidv.py
  compare_darkprop_iteration_flux.py

external/
  darkprop/
    modified DarkProp v0.3.0 source used for consistency checks
```


## Python setup

Use Python 3.10 or newer.

```text
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On Linux or macOS, activate the environment with:

```text
source .venv/bin/activate
```

## Quick checks

Validate the transformed iteration integrals:

```text
python codes/verify_iteration_integrals.py
```

Run a small smoke calculation:

```text
python codes/run_dirac_dm_isoscalar_vector_coupling.py --preset smoke --mchi-mev 100 --sigma-chin-cm2 1e-33
```

The default output location is `codes/output/`. This directory is ignored by
Git because it contains generated `.npz`, logs, and figures.

Plot cumulative detector spectra from a generated `.npz` file:

```text
python codes/plot_detector_cumulative_spectra.py codes/output/<run-output>.npz --x-axis speed
```

## DarkProp comparison workflow

The `external/darkprop/` tree is the modified DarkProp source copy that matches
the comparison scripts in `codes/`. Build products and Monte Carlo output
directories are excluded from version control.

Generate benchmark surface flux input files for DarkProp:

```text
python codes/generate_darkprop_surface_flux.py --mchi-mev 5000 --sigma-chin-cm2 1e-33
```

After running the corresponding DarkProp Monte Carlo, export and compare spectra
with:

```text
python codes/export_darkprop_mc_dphidv.py --darkprop <darkprop-output.hdf5>
python codes/compare_darkprop_iteration_flux.py --darkprop <darkprop-output.hdf5> --iteration <iteration-output.npz>
```

## Source relationship

This repository is not meant to be a full mirror of the manuscript working
directory. It is the code-release subset:

- `codes/*.py` contain the deterministic Boltzmann-iteration calculation and
  analysis helpers.
- `external/darkprop/` is copied from the modified DarkProp v0.3.0 source tree
  used for the Monte Carlo comparison.
- `codes/output/`, DarkProp build folders, and DarkProp example output folders
  are deliberately ignored.

## Citation

If you use these codes, cite the associated work:

Chuan-Yang Xing and Chen Xia, "The Boltzmann Equation Approach to the Dark
Matter Attenuation inside Earth".

## License

This code release is distributed under the MIT License. See `LICENSE`.
