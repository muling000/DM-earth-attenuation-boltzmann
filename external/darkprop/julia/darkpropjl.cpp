#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <jlcxx/const_array.hpp>
#include "darkprop/HomoEarth.hpp"
#include "darkprop/HomoElectronEarth.hpp"
#include "darkprop/PREMEarth.hpp"
#include "darkprop/Sun.hpp"
#include "darkprop/DMElectron.hpp"
#include "darkprop/SIDM.hpp"
#include "darkprop/SolarDM.hpp"
#include "darkprop/Simulation.hpp"

namespace dp = darkprop;
using Vector3d = dp::Vector3d<double>;
using RandomNumber = dp::RandomNumber<double>;
using Event = dp::Event<Vector3d, double>;
using Particle = dp::Particle<Vector3d, double>;
using ParticleElastic = dp::ParticleElastic<Vector3d, double>;
using Medium = dp::Medium<Vector3d, double>;
using MediumBall = dp::MediumBall<Vector3d, double>;
using Earth = dp::Earth<Vector3d, double>;
using HomoEarth = dp::HomoEarth<Vector3d, double>;
using HomoElectronEarth = dp::HomoElectronEarth<Vector3d, double>;
using PREMEarth = dp::PREMEarth<Vector3d, double>;
using Sun = dp::Sun<Vector3d, double>;
using DMElectron = dp::DMElectron<Vector3d, double>;
using SIDM = dp::SIDM<Vector3d, double>;
using SolarDM = dp::SolarDM<Vector3d, double>;

#ifdef __WITH_GENIE__
#include "darkprop/GENIEBDM.hpp"
#include "Framework/ParticleData/BaryonResUtils.h"
using GENIEBDM = dp::GENIEBDM<Vector3d, double>;
#endif

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
  mod.add_type<dp::Target>("Target")
    .constructor<const std::string&, int, int, double>()
    .method("get_Z", [](const dp::Target& t) { return t.Z; })
    .method("get_A", [](const dp::Target& t) { return t.A; })
    .method("get_m", [](const dp::Target& t) { return t.m; })
    .method("get_name", [](const dp::Target& t) { return t.name; });

  mod.add_type<Vector3d>("Vector3d", jlcxx::julia_type("AbstractFloat", "Base"))
    .constructor<double, double, double>()
    .constructor<const Vector3d&>()
    .method("isApprox", &Vector3d::isApprox)
    .method("norm", &Vector3d::norm)
    .method("normalized", &Vector3d::normalized)
    .method("dot", &Vector3d::dot)
    .method("cross", &Vector3d::cross);

  mod.set_override_module(jl_base_module);
  mod.method("+", [](const Vector3d& a, const Vector3d& b) { return a + b; });
  mod.method("-", [](const Vector3d& a, const Vector3d& b) { return a - b; });
  mod.method("-", [](const Vector3d& a) { return -a; });
  mod.method("==", [](const Vector3d& a, const Vector3d& b) { return a == b; });
  mod.method("vec", [](const Vector3d& a) {
      return std::vector<double> {a[0], a[1], a[2]};
    });
  mod.unset_override_module();

  mod.add_type<RandomNumber>("RandomNumber")
    .constructor<int>()
    .method("uniform_phi", &RandomNumber::uniform_phi)
    .method("uniform_costheta", &RandomNumber::uniform_costheta)
    .method("uniform_xi", &RandomNumber::uniform_xi)
    .method("uniform_ab", &RandomNumber::uniform_ab)
    .method("normal", &RandomNumber::normal);

  mod.add_type<Event>("Event")
    .constructor<double, double, Vector3d, Vector3d, double>()
    .method("get_t", [](const Event& e) { return e.t; })
    .method("get_T", [](const Event& e) { return e.T; })
    .method("get_r", [](const Event& e) { return e.r; })
    .method("get_p3", [](const Event& e) { return e.p3; })
    .method("get_weight", [](const Event& e) { return e.weight; });

  mod.add_type<Particle>("Particle")
    .method("pFromT", &Particle::pFromT)
    .method("TFromP3", &Particle::TFromP3)
    .method("updateEr!", &Particle::updateEr)
    .method("updateVwithP3T!", &Particle::updateVwithP3T)
    .method("updateP3TwithV!", &Particle::updateP3TwithV)
    .method("setR!", &Particle::setR)
    .method("setP!", &Particle::setP)
    .method("setV!", &Particle::setV)
    .method("totalCrossSection", &Particle::totalCrossSection)
    .method("scatter!", &Particle::scatter)
    .method("toEvent", &Particle::toEvent)
    .method("get_m", [](const Particle& p) { return p.m; })
    .method("get_sigma0", [](const Particle& p) { return p.sigma0; })
    .method("get_t", [](const Particle& p) { return p.t; })
    .method("get_T", [](const Particle& p) { return p.T; })
    .method("get_E", [](const Particle& p) { return p.T + p.m; })
    .method("get_r", [](const Particle& p) { return p.r; })
    .method("get_v", [](const Particle& p) { return p.v; })
    .method("get_p3", [](const Particle& p) { return p.p3; })
    .method("get_ep", [](const Particle& p) { return p.ep; })
    .method("get_er", [](const Particle& p) { return p.er; })
    .method("in_medium", [](const Particle& p) { return p.in_medium; })
    .method("in_medium!", [](Particle& p, bool in) { p.in_medium = in; });

  mod.add_type<ParticleElastic>("ParticleElastic", jlcxx::julia_base_type<Particle>())
    .method("scatteringCosTheta", &ParticleElastic::scatteringCosTheta)
    .method("scatter!", &ParticleElastic::scatter)
    .method("maximumRecoilT", &ParticleElastic::maximumRecoilT);

  mod.add_type<SIDM>("SIDM", jlcxx::julia_base_type<ParticleElastic>())
    .constructor<double, double, double>()
    // .constructor<double, double, double, Vector3d, Vector3d>()  // FIXME
    .method("sigmaEff", &SIDM::sigmaEff)
    .method("differentialCrossSection", &SIDM::differentialCrossSection)
    .method("totalCrossSection", &SIDM::totalCrossSection)
    .method("scatteringSampleTr", &SIDM::scatteringSampleTr);

  mod.add_type<DMElectron>("DMElectron", jlcxx::julia_base_type<ParticleElastic>())
    .constructor<double, double, double, double>()
    .method("FDMintegral", &DMElectron::FDMintegral)
    .method("inverseCDF", &DMElectron::inverseCDF)
    .method("totalCrossSection", &DMElectron::totalCrossSection)
    .method("scatteringSampleTr", &DMElectron::scatteringSampleTr);

  mod.add_type<SolarDM>("SolarDM", jlcxx::julia_base_type<Particle>())
    .constructor<double, double, double>()
    .method("setTempScale!", &SolarDM::setTempScale)
    .method("totalCrossSection", &SolarDM::totalCrossSection)
    .method("scatter!", &SolarDM::scatter);

  // Add the base class before GENIE
  mod.add_type<Medium>("Medium")
    .method("setCache!", &Medium::setCache)
    .method("getTargets", &Medium::getTargets)
    .method("propagate!", &Medium::propagate)
    .method("sampleTarget", &Medium::sampleTarget);

#ifdef __WITH_GENIE__
  mod.add_type<GENIEBDM>("GENIEBDM", jlcxx::julia_base_type<Particle>())
    .constructor<double, double, double, double, double>()
    .method("totalCrossSection",
      static_cast<double (GENIEBDM::*)(const dp::Target&) const>(
        &GENIEBDM::totalCrossSection))
    .method("totalCrossSection",
      static_cast<double (GENIEBDM::*)(const std::string&)>(&GENIEBDM::totalCrossSection))
    .method("scatter!", &GENIEBDM::scatter)
    .method("initGEVGDriver!",
      static_cast<void (GENIEBDM::*)(int, int, std::size_t)>(&GENIEBDM::initGEVGDriver))
    .method("initGEVGDriver!",
      static_cast<void (GENIEBDM::*)(const Medium&, std::size_t)>(&GENIEBDM::initGEVGDriver))
    .method("setgz!", static_cast<void (GENIEBDM::*)(void)>(&GENIEBDM::setgz))
    .method("setgz!", static_cast<void (GENIEBDM::*)(double)>(&GENIEBDM::setgz));

  mod.method("genie_init!", &dp::genie_init);
  mod.method("genie_set_coupling!", &dp::genie_set_coupling);
  mod.method("genie_add_dark_matter!",
      static_cast<void (*)(double, double, double)>(&dp::genie_add_dark_matter));
  mod.method("genie_add_dark_matter!",
      static_cast<void (*)(double, double)>(&dp::genie_add_dark_matter));
  mod.method("genie_set_xsec_table!", &dp::genie_set_xsec_table);
  mod.method("genie_set_event_generator_list!", &dp::genie_set_event_generator_list);
  mod.method("genie_reconfiguration!", &dp::genie_reconfiguration);
  mod.method("genie_reset_tune!", &dp::genie_reset_tune);
#endif

  mod.add_type<MediumBall>("MediumBall", jlcxx::julia_base_type<Medium>())
    .method("getRadius", &MediumBall::getRadius)
    .method("setRadius!", &MediumBall::setRadius)
    .method("getRFinal", &MediumBall::getRFinal)
    .method("setRFinal!", &MediumBall::setRFinal);

  mod.add_type<Earth>("Earth", jlcxx::julia_base_type<MediumBall>());

  mod.add_type<HomoEarth>("HomoEarth", jlcxx::julia_base_type<Earth>())
    .constructor<>()
    .method("setCache!", &HomoEarth::setCache)
    .method("meanFreePath", static_cast<double (HomoEarth::*)(const Particle&, std::size_t) const>(&HomoEarth::meanFreePath))
    .method("meanFreePath", static_cast<double (HomoEarth::*)(const Particle&) const>(&HomoEarth::meanFreePath))
    .method("sampleTarget", &HomoEarth::sampleTarget)
    .method("propagate!", &HomoEarth::propagate);

  mod.add_type<HomoElectronEarth>("HomoElectronEarth", jlcxx::julia_base_type<Earth>())
    .constructor<>()
    .method("setCache!", &HomoElectronEarth::setCache)
    .method("meanFreePath", &HomoElectronEarth::meanFreePath)
    .method("propagate!", &HomoElectronEarth::propagate)
    .method("sampleTarget", &HomoElectronEarth::sampleTarget);

  mod.add_type<PREMEarth>("PREMEarth", jlcxx::julia_base_type<Earth>())
    .constructor<bool, double, std::size_t>()
    .method("setCache!", &PREMEarth::setCache)
    .method("setPREM!", &PREMEarth::setPREM)
    .method("getLayer", &PREMEarth::getLayer)
    .method("getLayerRadius", &PREMEarth::getLayerRadius)
    .method("propagate!", &PREMEarth::propagate)
    .method("sampleTarget", &PREMEarth::sampleTarget)
    .method("toNextBoundary", &PREMEarth::toNextBoundary)
    .method("densityIntegral", &PREMEarth::densityIntegral)
    .method("density", &PREMEarth::density);

  mod.add_type<Sun>("Sun", jlcxx::julia_base_type<MediumBall>())
    .constructor<double, double>()
    .method("init!", &Sun::init)
    .method("densityIntegral", &Sun::densityIntegral)
    .method("densityIntegralNormalized", &Sun::densityIntegralNormalized)
    .method("densityIntegralInverseNormalized", &Sun::densityIntegralInverseNormalized)
    .method("propagate!", &Sun::propagate)
    .method("sampleTarget", &Sun::sampleTarget);

  mod.method("scatter!", &dp::scatter<Vector3d, double>, true);
  mod.method("propagate!", &dp::propagate<Vector3d, double>, true);
  mod.method("reduce_m", &dp::reduce_m<double>);

  // Simulation.hpp
  mod.method("simulate_cross_sphere", &dp::simulate_cross_sphere<Vector3d, double>);
  mod.method("simulate_cross_depth", &dp::simulate_cross_depth<Vector3d, double>);
  mod.method("simulate_track", &dp::simulate_track<Vector3d, double>);
  mod.method("analyse_track", &dp::analyse_track<Vector3d, double>);

#ifdef __WITH_GENIE__
  mod.method("E3_E0_p0_m3_m4_cosa", dp::E3_E0_p0_m3_m4_cosa<double>);
  mod.method("E3_E1_m1_m2_m3_m4_cosa", dp::E3_E1_m1_m2_m3_m4_cosa<double>);
  mod.method("dXSec_dQ2_DMCEL", dp::dXSec_dQ2_DMCEL<double>);
  mod.method("dXSec_dTchi_DMCEL", dp::dXSec_dTchi_DMCEL<double>);
  mod.method("dXSec_dcostheta_DMCEL", dp::dXSec_dcostheta_DMCEL<double>);
  mod.method("dXSec_dQ2_DMEL", dp::dXSec_dQ2_DMEL<double>);
  mod.method("dXSec_dTchi_DMEL", dp::dXSec_dTchi_DMEL<double>);
  mod.method("dXSec_dcostheta_DMEL", dp::dXSec_dcostheta_DMEL<double>);
  mod.method("dXSec_dQ2_DMRES", dp::dXSec_dQ2_DMRES<double>);
  mod.method("dXSec_dTchi_DMRES", dp::dXSec_dTchi_DMRES<double>);
  mod.method("dXSec_dcostheta_DMRES", dp::dXSec_dcostheta_DMRES<double>);
  mod.method("dXSec_dQ2_DMDIS", dp::dXSec_dQ2_DMDIS<double>);
  mod.method("dXSec_dTchi_DMDIS", dp::dXSec_dTchi_DMDIS<double>);
  mod.method("dXSec_dcostheta_DMDIS", dp::dXSec_dcostheta_DMDIS<double>);
  mod.method("total_xsec_DMCEL", dp::total_xsec_DMCEL<double>);
  mod.method("total_xsec_DMEL", dp::total_xsec_DMEL<double>);
  mod.method("total_xsec_DMRES", dp::total_xsec_DMRES<double>);
  mod.method("total_xsec_DMDIS", dp::total_xsec_DMDIS<double>);
  mod.method("limits_Echip_DMRES", dp::limits_Echip_DMRES<double>);
  mod.method("limits_Echip_DMDIS", dp::limits_Echip_DMDIS<double>);
  mod.method("limits_costheta_Echip", dp::limits_costheta_Echip<double>);
  mod.method("limits_W", dp::limits_W<double>);
  mod.method("limits_Q2", dp::limits_Q2<double>);
  mod.method("limits_Q2_W", dp::limits_Q2_W<double>);

  mod.add_bits<genie::EKinePhaseSpace>("KinePhaseSpace_t", jlcxx::julia_type("CppEnum"));
  mod.add_bits<genie::EResonance>("Resonance_t", jlcxx::julia_type("CppEnum"));
  mod.set_const("kNoResonance", genie::kNoResonance);
  mod.set_const("kP33_1232", genie::kP33_1232);
  mod.set_const("kS11_1535", genie::kS11_1535);
  mod.set_const("kD13_1520", genie::kD13_1520);
  mod.set_const("kS11_1650", genie::kS11_1650);
  mod.set_const("kD13_1700", genie::kD13_1700);
  mod.set_const("kD15_1675", genie::kD15_1675);
  mod.set_const("kS31_1620", genie::kS31_1620);
  mod.set_const("kD33_1700", genie::kD33_1700);
  mod.set_const("kP11_1440", genie::kP11_1440);
  mod.set_const("kP33_1600", genie::kP33_1600);
  mod.set_const("kP13_1720", genie::kP13_1720);
  mod.set_const("kF15_1680", genie::kF15_1680);
  mod.set_const("kP31_1910", genie::kP31_1910);
  mod.set_const("kP33_1920", genie::kP33_1920);
  mod.set_const("kF35_1905", genie::kF35_1905);
  mod.set_const("kF37_1950", genie::kF37_1950);
  mod.set_const("kP11_1710", genie::kP11_1710);
  mod.set_const("kF17_1970", genie::kF17_1970);

  mod.method("utilsres_AngularMom", genie::utils::res::AngularMom);
  mod.method("utilsres_Isospin", genie::utils::res::Isospin);
  mod.method("utilsres_Parity", genie::utils::res::Parity);
  mod.method("utilsres_AsString", [](genie::Resonance_t res) {
      return std::string(genie::utils::res::AsString(res)); });
  mod.method("utilsres_FromString", genie::utils::res::FromString);
  mod.method("utilsres_PdgCode", genie::utils::res::PdgCode);
  mod.method("utilsres_FromPdgCode", genie::utils::res::FromPdgCode);
  mod.method("utilsres_IsBaryonResonance", genie::utils::res::IsBaryonResonance);
  mod.method("utilsres_IsDelta", genie::utils::res::IsDelta);
  mod.method("utilsres_IsN", genie::utils::res::IsN);
  mod.method("utilsres_Mass", genie::utils::res::Mass);
  mod.method("utilsres_Width", genie::utils::res::Width);
  mod.method("utilsres_BWNorm", genie::utils::res::BWNorm);
  mod.method("utilsres_OrbitalAngularMom", genie::utils::res::OrbitalAngularMom);
  mod.method("utilsres_ResonanceIndex", genie::utils::res::ResonanceIndex);
  mod.method("utilsres_Cjsgn_plus", genie::utils::res::Cjsgn_plus);
  mod.method("utilsres_Dsgn", genie::utils::res::Dsgn);

  mod.add_type<genie::Algorithm>("Algorithm");

  auto algfactory = mod.add_type<genie::AlgFactory>("AlgFactory");

  mod.method("InstanceAlgFactory", &genie::AlgFactory::Instance);

  mod.add_type<genie::Interaction>("Interaction")
    .constructor<>();

  mod.method("Interaction_DMRES", [](int tgt, int nuc, int probe, double E) {
          return genie::Interaction::DMRES(tgt, nuc, probe, E);
        });

  mod.add_type<genie::XSecAlgorithmI>("XSecAlgorithmI",
      jlcxx::julia_base_type<genie::Algorithm>())
    .method("XSec", &genie::XSecAlgorithmI::XSec)
    .method("Integral", &genie::XSecAlgorithmI::Integral)
    .method("ValidProcess", &genie::XSecAlgorithmI::ValidProcess)
    .method("ValidKinematics", &genie::XSecAlgorithmI::ValidKinematics);


  mod.add_type<genie::DMRESPXSec>("DMRESPXSec",
    jlcxx::julia_base_type<genie::XSecAlgorithmI>())
    .constructor<>()
    .method("XSec", &genie::DMRESPXSec::XSec)
    .method("Integral", &genie::DMRESPXSec::Integral)
    .method("Configure!", static_cast<void (genie::DMRESPXSec::*)(std::string)>(
          &genie::DMRESPXSec::Configure))
    .method("F1u", &genie::DMRESPXSec::F1u)
    .method("F1d", &genie::DMRESPXSec::F1d)
    .method("F2u", &genie::DMRESPXSec::F2u)
    .method("F2d", &genie::DMRESPXSec::F2d)
    .method("FAu", &genie::DMRESPXSec::FAu)
    .method("FAd", &genie::DMRESPXSec::FAd)
    .method("FPu", &genie::DMRESPXSec::FPu)
    .method("FPd", &genie::DMRESPXSec::FPd)
    .method("C3Au", &genie::DMRESPXSec::C3Au)
    .method("C3Ad", &genie::DMRESPXSec::C3Ad)
    .method("C4Au", &genie::DMRESPXSec::C4Au)
    .method("C4Ad", &genie::DMRESPXSec::C4Ad)
    .method("C5Au", &genie::DMRESPXSec::C5Au)
    .method("C5Ad", &genie::DMRESPXSec::C5Ad)
    .method("C6Au", &genie::DMRESPXSec::C6Au)
    .method("C6Ad", &genie::DMRESPXSec::C6Ad)
    .method("C3Vu", &genie::DMRESPXSec::C3Vu)
    .method("C3Vd", &genie::DMRESPXSec::C3Vd)
    .method("C4Vu", &genie::DMRESPXSec::C4Vu)
    .method("C4Vd", &genie::DMRESPXSec::C4Vd)
    .method("C5Vu", &genie::DMRESPXSec::C5Vu)
    .method("C5Vd", &genie::DMRESPXSec::C5Vd)
    .method("C6Vu", &genie::DMRESPXSec::C6Vu)
    .method("C6Vd", &genie::DMRESPXSec::C6Vd)
    .method("F1", &genie::DMRESPXSec::F1)
    .method("F2", &genie::DMRESPXSec::F2)
    .method("C3V", &genie::DMRESPXSec::C3V)
    .method("C4V", &genie::DMRESPXSec::C4V)
    .method("C5V", &genie::DMRESPXSec::C5V)
    .method("C6V", &genie::DMRESPXSec::C6V)
    .method("FA", &genie::DMRESPXSec::FA)
    .method("FP", &genie::DMRESPXSec::FP)
    .method("C3A", &genie::DMRESPXSec::C3A)
    .method("C4A", &genie::DMRESPXSec::C4A)
    .method("C5A", &genie::DMRESPXSec::C5A)
    .method("C6A", &genie::DMRESPXSec::C6A)
    .method("Wcut", &genie::DMRESPXSec::Wcut)
    .method("Q2cut", &genie::DMRESPXSec::Q2cut);

  algfactory.method("GetAlgorithm", static_cast<
        const genie::Algorithm* (genie::AlgFactory::*)(string, string)>(
          &genie::AlgFactory::GetAlgorithm))
    .method("get_DMRESPXSec",
      [](const genie::AlgFactory *alg, const std::string& name, const std::string& param) {
        return dynamic_cast<const genie::DMRESPXSec *>(
          genie::AlgFactory::Instance()->GetAlgorithm(name, param));
      });
#endif
}

namespace jlcxx
{
  template<> struct SuperType<ParticleElastic> { using type = Particle; };
  template<> struct SuperType<SIDM> { using type = ParticleElastic; };
  template<> struct SuperType<DMElectron> { using type = ParticleElastic; };
  template<> struct SuperType<MediumBall> { using type = Medium; };
  template<> struct SuperType<Earth> { using type = MediumBall; };
  template<> struct SuperType<HomoEarth> { using type = Earth; };
  template<> struct SuperType<HomoElectronEarth> { using type = Earth; };
  template<> struct SuperType<PREMEarth> { using type = Earth; };
  template<> struct SuperType<Sun> { using type = MediumBall; };
#ifdef __WITH_GENIE__
  template<> struct SuperType<GENIEBDM> { using type = Particle; };
  template<> struct SuperType<genie::XSecAlgorithmI> { using type = genie::Algorithm; };
  template<> struct SuperType<genie::DMRESPXSec> { using type = genie::XSecAlgorithmI; };
#endif
}
