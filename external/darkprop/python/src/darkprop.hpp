#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "darkprop/core.hpp"
namespace py = pybind11;
namespace dp = darkprop;
using Vector3d = dp::Vector3d<double>;
using RandomNumber = dp::RandomNumber<double>;
using Particle = dp::Particle<Vector3d, double>;
using Earth = dp::Earth<Vector3d, double>;
using namespace pybind11::literals;
