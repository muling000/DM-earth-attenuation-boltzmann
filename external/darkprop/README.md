# DarkProp

Monte Carlo simulation code for the propagation of dark matter particles in a medium.


# Installation

DarkProp provides an executable application `darkprop` based on the darkprop C++ library.
The `darkprop` application reads a [TOML](https://toml.io/) configuration file as input and
stores outputs in [HDF5](https://www.hdfgroup.org/) format. The numerical calculations
require [GSL](https://www.gnu.org/software/gsl/) library, and the simulation is
parallelized using [MPI](https://www.mpi-forum.org/).
The [spdlog](https://github.com/gabime/spdlog) library is used for logging.

Above mentioned dependencies and the build tools (`cmake>=3.20`, `g++`/`clang` with `C++17`
support) are required to build the application. They are available in the repositories of
mainstream Linux distributions and from [homebrew](https://brew.sh/) for macOS, for example

- Ubuntu 22.04+

        $ sudo apt install g++ cmake libgsl-dev libhdf5-mpi-dev libspdlog-dev

- Fedora 33+

        $ sudo dnf install g++ cmake mpich-autoload gsl-devel hdf5-mpich-devel spdlog-devel
        $ source /etc/profile

- macOS (arm64)

        $ brew install cmake gsl hdf5-mpi spdlog

The source code of the `yuc::config` tool for parsing TOML configurations and the
[indicators](https://github.com/p-ranav/indicators) library for showing progress bar have
been included in the source.

After installing the dependencies, `darkprop` can be build with `cmake` as follows.

Get into the source

    $ tar -xzvf darkprop-vX.Y.Z.tar.gz
    $ cd darkprop-vX.Y.Z
    $ mkdir build; cd build

Configure

    $ cmake ..

Build and install

    $ cmake --build . --target install

The default location of installation may require root privilege, in which case the last
step should be run with `sudo`. Custom location can be specified by setting
`-DCMAKE_INSTALL_PREFIX` at the configuration step, e.g.

    $ cmake -DCMAKE_INSTALL_PREFIX=/path/to/a/location ..

In this case, you may need to add the `/path/to/a/location/bin` path to the `PATH`
environment variable to make the `darkprop` command available.

# Uninstallation

In the `build` directory where the `cmake` command was executed, type

    $ cmake --build . --target uninstall

# Usage

The `darkprop` application can be launched as follows,

    $ mpiexec -n 4 darkprop config.toml

where 4 is the number of CPU cores for parallel calculation. It can be replaced with any
other number that fits your machine. `config.toml` is the configuration file. A complete
example can be found in the `app` directory. The explanation of the parameters is
documented at [https://darkprop.site](https://darkprop.site).

After the simulation, the crossing events at underground surfaces will be stored.
The reconstruction of the underground DM flux is described in the paper [^1] and is
currently not integrated in the `darkprop` application.

An example of the Earth attenuation of isotropic CRDM flux is provided in the
`examples/crdm-isotropic` directory, where the statistical analysis is implemented with
Python.


## C++ library

The core of DarkProp is the darkprop C++ library, which can be used by including the header

    #include "darkprop/core.hpp"
    #include "darkprop/IO.hpp"  // optional hdf5 io

See the example `examples/trajectory-csv` for basic usage.


## Python

We also provide a Python package which is a wrapper of the darkprop C++ library.
It can be installed via `pip`

    $ pip install darkprop

See `examples/python/basic-example.ipynb` for more information.


# Citation

If you use this code in your published work, please cite our paper [^1]
([`CITATION.bib`](CITATION.bib)).

[^1]: C. Xia, Y. H. Xu, Y. F. Zhou, "Production and attenuation of cosmic-rays
      boosted dark matter". JCAP, 2022, 02(02): 028, arXiv: 2111.05559
