#include "darkprop.hpp"
#include "darkprop/Injector.hpp"

using Injector = dp::Injector<Vector3d, double>;
using HaloInjector = dp::HaloInjector<Vector3d, double>;
using FluxInjector = dp::FluxInjector<Vector3d, double>;
using IntensityInjector = dp::IntensityInjector<Vector3d, double>;
using SampleInjector = dp::SampleInjector<Vector3d, double>;
using SourceInjector = dp::SourceInjector<Vector3d, double>;

template<typename InjectorBase=Injector> class PyInjector : public InjectorBase
{
public:
  using InjectorBase::InjectorBase;
  double next(Particle& p) override {
    PYBIND11_OVERRIDE_PURE(double, InjectorBase, next, p);
  }
};

void bind_Injector(py::module_ &m)
{
  py::class_<Injector, PyInjector<>>(m, "Injector")
    .def("next", &Injector::next)
    .def("__call__", &Injector::operator())
    .def("setSeed", &Injector::setSeed)
    .def("sett0", &Injector::sett0);

  py::class_<HaloInjector, Injector>(m, "HaloInjector")
    .def(py::init<double, double, Vector3d, double, double, double, double, int>(),
        "vmin"_a=0,
        "vmax"_a=-1,
        "vearth"_a=Vector3d{0, 0, 240 * dp::units::km / dp::units::sec},
        "vesc"_a=dp::constants::vesc,
        "v0"_a=dp::constants::v0,
        "radius"_a=dp::constants::rEarth,
        "t0"_a=0,
        "seed"_a=-1)
    .def("next", &HaloInjector::next)
    .def("set_norm", &HaloInjector::set_norm);

  py::class_<FluxInjector, Injector>(m, "FluxInjector")
    .def(py::init<int, double>(), "seed"_a=-1, "radius"_a=dp::constants::rEarth)
    .def(py::init<const std::vector<double>&, const std::vector<double>&, double, int>(),
        "T"_a, "flux"_a, "radius"_a=dp::constants::rEarth, "seed"_a=-1)
    .def("init", &FluxInjector::init)
    .def("next", &FluxInjector::next);

  py::class_<IntensityInjector, Injector>(m, "IntensityInjector")
    .def(py::init<int, double>(), "seed"_a=-1, "radius"_a=dp::constants::rEarth)
    .def(py::init<const std::vector<double>&, const std::vector<double>&,
                  const std::vector<double>&, const std::vector<double>&,
                  int, double>(),
        "b"_a, "l"_a, "T"_a, "I"_a, "seed"_a=-1, "radius"_a=dp::constants::rEarth)
    .def("init", &IntensityInjector::init)
    .def("setSigma", &IntensityInjector::setSigma)
    .def("setBaseSigma", &IntensityInjector::setBaseSigma)
    .def("getNorm", &IntensityInjector::getNorm)
    .def("intensity", &IntensityInjector::intensity)
    .def("next", &IntensityInjector::next);

  py::class_<SampleInjector, Injector>(m, "SampleInjector")
    .def(py::init<int, double>(), "seed"_a=-1, "radius"_a=dp::constants::rEarth)
    .def(py::init<const std::string&, unsigned, int, double>(), "filename"_a,
        "flag"_a=H5F_ACC_RDONLY, "seed"_a=-1, "radius"_a=dp::constants::rEarth)
    .def("next", &SampleInjector::next);

  py::class_<SourceInjector, Injector>(m, "SourceInjector")
    .def(py::init<int>(), "seed"_a=-1)
    .def(py::init<const std::vector<double>&, const std::vector<double>&,
                  const std::vector<double>&, const std::vector<double>&,
                  const std::string&, int>(),
        "xi_T"_a, "inv_cdf_T"_a, "xi_r"_a, "inv_cdf_r"_a, "method"_a="linear", "seed"_a=-1)
    .def(py::init<const std::string&, const std::string&, const std::string&, int>(),
        "r_file"_a, "T_file"_a, "method"_a="linear", "seed"_a=-1)
    .def("interpolate_r", &SourceInjector::interpolate_r)
    .def("interpolate_T", &SourceInjector::interpolate_T)
    .def("next", &SourceInjector::next);
}
