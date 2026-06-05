#include "darkprop.hpp"
#include "darkprop/ParticleElastic.hpp"
#include "darkprop/SIDM.hpp"
#include "darkprop/SIHelmDM.hpp"
#include "darkprop/DMElectron.hpp"
#include "darkprop/SolarDM.hpp"

using ParticleElastic = dp::ParticleElastic<Vector3d, double>;
using SIDM = dp::SIDM<Vector3d, double>;
using SIHelmDM = dp::SIHelmDM<Vector3d, double>;
using DMElectron = dp::DMElectron<Vector3d, double>;
using SolarDM = dp::SolarDM<Vector3d, double>;
using Medium = dp::Medium<Vector3d, double>;

#ifdef __WITH_GENIE__
#include "darkprop/GENIEBDM.hpp"
using GENIEBDM = dp::GENIEBDM<Vector3d, double>;
#endif

template<typename ParticleBase=Particle> class PyParticle : public ParticleBase
{
public:
  using ParticleBase::ParticleBase;
  double scatter(const dp::Target& t, dp::RandomNumber<double>& rn) override {
    PYBIND11_OVERRIDE_PURE(double, ParticleBase, scatter, t, rn);
  }
  double totalCrossSection(const dp::Target& t) const override {
    PYBIND11_OVERRIDE_PURE(double, ParticleBase, totalCrossSection, t);
  }
};

template<typename ParticleElasticBase=ParticleElastic>
class PyParticleElastic : public PyParticle<ParticleElasticBase>
{
public:
  using PyParticle<ParticleElasticBase>::PyParticle;
  double scatteringSampleTr(const dp::Target& t,
                            dp::RandomNumber<double>& rn) override {
    PYBIND11_OVERRIDE_PURE(double, ParticleElasticBase, scatteringSampleTr, t, rn);
  }
  double scatteringCosTheta(const dp::Target& t, double Tr) override {
    PYBIND11_OVERRIDE(double, ParticleElasticBase, scatteringCosTheta, t, Tr);
  }
  double maximumRecoilT(const dp::Target& t) const override {
    PYBIND11_OVERRIDE(double, ParticleElasticBase, maximumRecoilT, t);
  }
  double scatter(const dp::Target& t, dp::RandomNumber<double>& rn) override {
    PYBIND11_OVERRIDE(double, ParticleElasticBase, scatter, t, rn);
  }
};

void bind_Particle(py::module_ &m)
{
  py::class_<Particle, PyParticle<>>(m, "Particle")
    .def_readwrite("m", &Particle::m)
    .def_readwrite("sigma0", &Particle::sigma0)
    .def_readwrite("t", &Particle::t)
    .def_readwrite("T", &Particle::T)
    .def_readwrite("r", &Particle::r)
    .def_readwrite("v", &Particle::v)
    .def_readwrite("p3", &Particle::p3)
    .def_readwrite("ep", &Particle::ep)
    .def_readwrite("er", &Particle::er)
    .def_readwrite("in_medium", &Particle::in_medium)
    .def("pFromT", &Particle::pFromT)
    .def("TFromP3", &Particle::TFromP3)
    .def("updateEr", &Particle::updateEr)
    .def("updateEp", &Particle::updateEp)
    .def("updateVwithP3T", &Particle::updateVwithP3T)
    .def("updateP3TwithV", &Particle::updateP3TwithV)
    .def("setR", &Particle::setR)
    .def("setP", &Particle::setP)
    .def("setV", &Particle::setV)
    .def("toEvent", &Particle::toEvent, "weight"_a=1.0, "kmsec"_a=true)
    .def("totalCrossSection", &Particle::totalCrossSection)
    .def("scatter", &Particle::scatter);

  py::class_<ParticleElastic, Particle, PyParticleElastic<>>(m, "ParticleElastic")
    .def("scatteringSampleTr", &ParticleElastic::scatteringSampleTr)
    .def("scatteringCosTheta", &ParticleElastic::scatteringCosTheta)
    .def("maximumRecoilT", &ParticleElastic::maximumRecoilT)
    .def("scatter", &ParticleElastic::scatter);

  py::class_<SIDM, ParticleElastic, PyParticleElastic<SIDM>>(m, "SIDM")
    .def(py::init<double, double, double>(),
         "sigma0"_a=1e-30*dp::units::cm2,
         "m"_a=1.0*dp::units::GeV,
         "t"_a=0.0)
    .def("sigmaEff", &SIDM::sigmaEff)
    .def("differentialCrossSection", &SIDM::differentialCrossSection)
    .def("totalCrossSection", &SIDM::totalCrossSection)
    .def("scatteringSampleTr", &SIDM::scatteringSampleTr);

  py::class_<SIHelmDM, SIDM>(m, "SIHelmDM")
    .def(py::init<double, double, double>(),
         "sigma0"_a=1e-30*dp::units::cm2,
         "m"_a=1.0*dp::units::GeV,
         "t"_a=0.0)
    .def("differentialCrossSection", &SIHelmDM::differentialCrossSection)
    .def("initCrossSectionCDF", &SIHelmDM::initCrossSectionCDF)
    .def("ff2Integral", &SIHelmDM::ff2Integral)
    .def("inverseCDF", &SIHelmDM::inverseCDF)
    .def("scatteringCosTheta", &SIHelmDM::scatteringCosTheta)
    .def("totalCrossSection", &SIHelmDM::totalCrossSection)
    .def("scatteringSampleTr", &SIHelmDM::scatteringSampleTr);

  py::class_<DMElectron, ParticleElastic>(m, "DMElectron")
    .def(py::init<double, double, double, double>(),
         "sigma0"_a=1e-30*dp::units::cm2,
         "m"_a=1.0*dp::units::GeV,
         "mmed"_a=1.0*dp::units::GeV,
         "t"_a=0.0)
    .def_readwrite("mediator_limit", &DMElectron::mediator_limit)
    .def_readwrite("abstol", &DMElectron::abstol)
    .def_readwrite("reltol", &DMElectron::reltol)
    .def_readwrite("maxit", &DMElectron::maxit)
    .def("FDM", &DMElectron::FDM)
    .def("FDMintegral", &DMElectron::FDMintegral)
    .def("inverseCDF", &DMElectron::inverseCDF)
    .def("totalCrossSection", &DMElectron::totalCrossSection)
    .def("scatteringSampleTr", &DMElectron::scatteringSampleTr);

  py::class_<SolarDM, Particle>(m, "SolarDM")
    .def(py::init<double, double, double>(),
         "sigma0"_a=1e-30*dp::units::cm2,
         "m"_a=1.0*dp::units::keV,
         "t"_a=0.0)
    .def("setTempScale", &SolarDM::setTempScale)
    .def("getTempScale", &SolarDM::getTempScale)
    .def("totalCrossSection", &SolarDM::totalCrossSection)
    .def("scatter", &SolarDM::scatter);

#ifdef __WITH_GENIE__
  m.def("genie_init", &dp::genie_init);
  m.def("genie_set_coupling", &dp::genie_set_coupling);
  m.def("genie_add_dark_matter",
      py::overload_cast<double, double, double>(&dp::genie_add_dark_matter));
  m.def("genie_add_dark_matter",
      py::overload_cast<double, double>(&dp::genie_add_dark_matter));
  m.def("genie_set_xsec_table", &dp::genie_set_xsec_table);
  m.def("genie_set_event_generator_list", &dp::genie_set_event_generator_list);
  m.def("genie_reconfiguration", &dp::genie_reconfiguration);

  py::class_<GENIEBDM, Particle>(m, "GENIEBDM")
    .def(py::init<double, double, double, double, double>(),
         "sigma0"_a=1e-30*dp::units::cm2,
         "t"_a=0.0,
         "Tmin"_a=1.1*dp::units::GeV,
         "Tmax"_a=10.0*dp::units::GeV,
         "gz_org"_a=1)
    .def("totalCrossSection", py::overload_cast<const dp::Target&>(
         &GENIEBDM::totalCrossSection, py::const_))
    .def("totalCrossSection", py::overload_cast<const std::string&>(
         &GENIEBDM::totalCrossSection))
    .def("scatter", &GENIEBDM::scatter)
    .def("initGEVGDriver", py::overload_cast<int, int, std::size_t>(&GENIEBDM::initGEVGDriver))
    .def("initGEVGDriver", py::overload_cast<const Medium&, std::size_t>(&GENIEBDM::initGEVGDriver));
#endif
}
