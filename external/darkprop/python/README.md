[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://en.wikipedia.org/wiki/MIT_License)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Python wrapper of DarkProp, a Monte Carlo simulation code for the propagation of dark
matter particles in a medium.

# Installation


## Install using pip

```bash
$ pip install --user darkprop
```
---
**Tip**

If there is a network problem during `pip install`, you can set up a mirror server.
For example, using [Tsinghua mirror](https://mirrors.tuna.tsinghua.edu.cn/help/pypi/),
add `-i` option to `pip install` command
```bash
$ pip install -i https://pypi.tuna.tsinghua.edu.cn/simple darkprop
```

---


## Install from the source

1. Download the source code from [DarkProp's homepage](http://yfzhou.itp.ac.cn/darkprop).

2. Install dependencies. A C++ compiler supporting `C++17` standard is required.
   [GSL](https://www.gnu.org/software/gsl/), [HDF5](https://www.hdfgroup.org/), and
   [spdlog](https://github.com/gabime/spdlog/) are needed to build the Python wrapper.
   They can be installed, for example

   * Ubuntu >= 22.04

       ```bash
       $ sudo apt install g++ libgsl-dev libhdf5-dev libspdlog-dev
       ```

   * Fedora >= 33

       ```bash
       $ sudo dnf install g++ python3-devel gsl-devel hdf5-devel spdlog-devel
       ```

   Then prepare a python virtual environment and activate it.
   ```bash
   $ python3 -m venv venv
   $ source venv/bin/activate
   ```

3. Build and install into the virtual environment.
   In the root directory of the source code, execute
   ```bash
   $ python3 -m pip install .
   ```


# Example

Install [Jupyter notebook](https://jupyter.org/) and [Matplotlib](https://matplotlib.org/)
to run the example.

```bash
$ jupyter-notebook example/python/basic-example.ipynb
```


# Uninstall

```bash
$ pip uninstall darkprop
```


# Citation

If you use darkprop in your publications, please cite the paper

- C. Xia, Y.-H. Xu, and Y.-F. Zhou, Production and attenuation of cosmic-ray
  boosted dark matter,
  [JCAP, 2022, 02(02): 028](https://doi.org/10.1088/1475-7516/2022/02/028),
  [arXiv: 2111.05559](https://arxiv.org/abs/2111.05559)
