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
  plot_damascus_boltzmann_dphidv.py

external/
  darkprop/
    modified DarkProp v0.3.0 source used for consistency checks
  damascus/
    modified DaMaSCUS source and Boltzmann benchmark configuration
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

## DaMaSCUS benchmark workflow

&emsp;&emsp;The `external/damascus/` tree contains the modified DaMaSCUS source used for the homogeneous-Earth benchmark. The benchmark-specific changes and parameters are summarized in `external/damascus/BOLTZMANN_BENCHMARK.md`.

Build and run DaMaSCUS from its source directory:

```text
cd external/damascus
make
cd bin
mpirun -n 4 ./DaMaSCUS-Simulator boltzmann_benchmark.cfg
mpirun -n 4 ./DaMaSCUS-Analyzer boltzmann_benchmark_5GeV_5e-32_vcut10_3m
```

&emsp;&emsp;Building DaMaSCUS requires an MPI C++ compiler, Eigen, libconfig++, and `pkg-config`. If Eigen is not discoverable through `pkg-config`, pass its compiler flag explicitly, for example `make EIGEN_CFLAGS=-I/path/to/eigen3`.

&emsp;&emsp;The full benchmark requests three million detector samples and is not a quick smoke test. Its generated `data/` and `results/` files are ignored by Git.

Export the fine-bin DarkProp table used by the three-way comparison:

```text
python codes/export_darkprop_mc_dphidv.py --darkprop <run1.hdf5> <run2.hdf5> --nbins 1000 --vmin 3.3356409519815205e-6 --output codes/output/darkprop_mc_dphidv_vmin1kms_fine1000.dat
```

After generating the deterministic Boltzmann result, create the three-way comparison with:

```text
python codes/plot_damascus_boltzmann_dphidv.py
```

&emsp;&emsp;Use `--damascus-data`, `--boltzmann`, and `--darkprop` to select non-default output locations. Generated plots are written to `codes/output/`.

## Source relationship

This repository is not meant to be a full mirror of the manuscript working
directory. It is the code-release subset:

- `codes/*.py` contain the deterministic Boltzmann-iteration calculation and
  analysis helpers.
- `external/darkprop/` is copied from the modified DarkProp v0.3.0 source tree
  used for the Monte Carlo comparison.
- `external/damascus/` is copied from the modified DaMaSCUS source tree used for
  the homogeneous-Earth Monte Carlo benchmark.
- `codes/output/` and Monte Carlo build and output folders are deliberately
  ignored.

## Citation

If you use these codes, cite the associated work ([`CITATION.bib`](CITATION.bib)):

Chuan-Yang Xing and Chen Xia, "Dark Matter Attenuation inside the Earth: A Boltzmann Equation Approach", arXiv:2606.16204.

```bibtex
@article{Xing:2026bbq,
    author = "Xing, Chuan-Yang and Xia, Chen",
    title = "{Dark Matter Attenuation inside the Earth: A Boltzmann Equation Approach}",
    eprint = "2606.16204",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    month = "6",
    year = "2026"
}
```

## License

&emsp;&emsp;This code release is distributed under the MIT License. See `LICENSE`. Bundled DarkProp and DaMaSCUS sources retain the license files included in their respective directories.
