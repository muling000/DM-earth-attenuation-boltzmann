#include "darkprop.hpp"
#include "darkprop/HomoEarth.hpp"
#include "darkprop/HomoElectronEarth.hpp"
#include "darkprop/PREMEarth.hpp"
#include "darkprop/Sun.hpp"

using Medium = dp::Medium<Vector3d, double>;
using MediumBall = dp::MediumBall<Vector3d, double>;
using HomoEarth = dp::HomoEarth<Vector3d, double>;
using HomoElectronEarth = dp::HomoElectronEarth<Vector3d, double>;
using PREMEarth = dp::PREMEarth<Vector3d, double>;
using Sun = dp::Sun<Vector3d, double>;

template<typename MediumBase=Medium> class PyMedium : public MediumBase
{
public:
  using MediumBase::MediumBase;
  void setCache(const Particle& t) override {
    PYBIND11_OVERRIDE(void, MediumBase, setCache, t);
  }
  const std::vector<dp::Target>& getTargets() const override {
    PYBIND11_OVERRIDE(const std::vector<dp::Target>&, MediumBase, getTargets,);
  }
  double propagate(Particle &p, RandomNumber& rn) override {
    PYBIND11_OVERRIDE_PURE(double, MediumBase, propagate, p, rn);
  }
  dp::Target sampleTarget(const Particle &p, RandomNumber &rn) const override {
    PYBIND11_OVERRIDE_PURE(dp::Target, MediumBase, sampleTarget, p, rn);
  }
};


void bind_Earth(py::module_ &m)
{
  py::class_<Medium, PyMedium<>>(m, "Medium")
    .def("setCache", &Medium::setCache)
    .def("getTargets", &Medium::getTargets)
    .def("propagate", &Medium::propagate)
    .def("sampleTarget", &Medium::sampleTarget);
  py::class_<MediumBall, Medium>(m, "MediumBall")
    .def("getRadius", &MediumBall::getRadius)
    .def("setRadius", &MediumBall::setRadius)
    .def("getRFinal", &MediumBall::getRFinal)
    .def("setRFinal", &MediumBall::setRFinal);
  py::class_<Earth, MediumBall>(m, "Earth");

  py::class_<HomoEarth, Earth>(m, "HomoEarth")
    .def(py::init<>())
    .def("setCache", &HomoEarth::setCache)
    .def("meanFreePath", py::overload_cast<const Particle&, std::size_t>(&HomoEarth::meanFreePath, py::const_))
    .def("meanFreePath", py::overload_cast<const Particle&>(&HomoEarth::meanFreePath, py::const_))
    .def("propagate", &HomoEarth::propagate)
    .def("sampleTarget", &HomoEarth::sampleTarget);

  py::class_<HomoElectronEarth, Earth>(m, "HomoElectronEarth")
    .def(py::init<>())
    .def("setCache", &HomoElectronEarth::setCache)
    .def("meanFreePath", &HomoElectronEarth::meanFreePath)
    .def("propagate", &HomoElectronEarth::propagate)
    .def("sampleTarget", &HomoElectronEarth::sampleTarget);

  py::class_<PREMEarth, Earth>(m, "PREMEarth")
    .def(py::init<bool, double, std::size_t>(),
        "energy_dependent"_a=false,
        "accuracy"_a=1e-12,
        "maxit"_a=100)
    .def("getLayer", &PREMEarth::getLayer)
    .def("getLayerRadius", &PREMEarth::getLayerRadius)
    .def("setPREM", &PREMEarth::setPREM)
    .def("setCache", &PREMEarth::setCache)
    .def("propagate", &PREMEarth::propagate)
    .def("sampleTarget", &PREMEarth::sampleTarget)
    .def("toNextBoundary", &PREMEarth::toNextBoundary)
    .def("densityIntegral", &PREMEarth::densityIntegral)
    .def("density", &PREMEarth::density);

  py::class_<Sun, MediumBall>(m, "Sun")
    .def(py::init<double, double>(), "r"_a=dp::constants::rSun, "rf"_a=dp::units::cm)
    .def(py::init<const std::string&, double, double>(), "filename"_a, "r"_a=dp::constants::rSun, "rf"_a=dp::units::cm)
    .def("init", &Sun::init)
    .def("densityIntegral", &Sun::densityIntegral)
    .def("densityIntegralNormalized", &Sun::densityIntegralNormalized)
    .def("densityIntegralInverseNormalized", &Sun::densityIntegralInverseNormalized)
    .def("propagate", &Sun::propagate)
    .def("sampleTarget", &Sun::sampleTarget);
}
