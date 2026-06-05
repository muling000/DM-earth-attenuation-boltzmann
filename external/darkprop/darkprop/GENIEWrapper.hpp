/**
 * @file
 * Tool functions for using the GENIE code.
 */
#pragma once
#include <Framework/Utils/AppInit.h>
#include <Framework/Utils/RunOpt.h>
#include <Framework/Utils/KineUtils.h>
#include <Framework/Algorithm/AlgConfigPool.h>
#include <Framework/ParticleData/PDGCodes.h>
#include <Framework/ParticleData/PDGLibrary.h>
#include <Framework/ParticleData/PDGUtils.h>
#include <Framework/GHEP/GHepRecord.h>
#include <Framework/Conventions/Units.h>
#include <Physics/BoostedDarkMatter/XSection/AhrensDMELPXSec.h>
#include <Physics/BoostedDarkMatter/XSection/QPMDMDISPXSec.h>
#include <Physics/BoostedDarkMatter/XSection/DMRESPXSec.h>
#include <Physics/BoostedDarkMatter/XSection/DMCELPXSec.h>
#include <Physics/BoostedDarkMatter/XSection/InelXSecUtils.h>
#include "darkprop/Base.hpp"

namespace darkprop {

inline void genie_reset_tune(const std::string& tune_name="GDM23_00a_01_000")
{
    genie::RunOpt::Instance()->SetTuneName(tune_name);
    genie::RunOpt::Instance()->BuildTune();
    genie::AlgConfigPool::Instance()->ForceDelete();
    genie::AlgFactory::Instance()->ForceDelete();
}

inline void genie_init(
    const std::string& tune_name="GDM23_00a_01_000",
    const std::string& mesg_thre="Messenger_whisper.xml",
    const std::string& xsec_file="",
    const std::string& evgen_list="DMHAD-NOFSI",
    int seed=0,
    const std::string& cache_file="")
{
    genie::RunOpt * runopt = genie::RunOpt::Instance();
    if (!runopt->Tune()) {
        runopt->SetTuneName(tune_name);
        runopt->BuildTune();
    } else {
        genie_reset_tune(tune_name);
    }
    runopt->EnableBareXSecPreCalc(false);
    runopt->SetEventGeneratorList(evgen_list);
    genie::utils::app_init::MesgThresholds(mesg_thre);
    genie::utils::app_init::CacheFile(cache_file);
    genie::utils::app_init::XSecTable(xsec_file, false);
    genie::utils::app_init::RandGen(seed);
    genie::GHepRecord::SetPrintLevel(runopt->EventRecordPrintLevel());
}

inline void genie_set_coupling(double gz)
{
    genie::Registry * r = genie::AlgConfigPool::Instance()->CommonList(
            "Param", "BoostedDarkMatter");
    r->UnLock();
    r->Set("ZpCoupling", gz);
    r->Lock();
    genie::AlgFactory::Instance()->ForceReconfiguration();
}

inline void genie_add_dark_matter(double xsec0, double mchi, double mmed)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    pdglib->ReloadDBase();
    pdglib->AddDarkMatter(mchi, mmed / mchi);
    double gz = mmed * std::pow(M_PI * xsec0
        / (9.0 * std::pow(reduce_m(mchi, constants::mp), 2.0)), 1.0/4.0);
    genie_set_coupling(gz);
}

inline void genie_add_dark_matter(double mchi, double mmed)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    pdglib->ReloadDBase();
    pdglib->AddDarkMatter(mchi, mmed / mchi);
    // Make sure that mediator mass is updated
    genie::AlgFactory::Instance()->ForceReconfiguration();
}

inline void genie_set_xsec_table(const std::string& t_xsec_file) {
    genie::utils::app_init::XSecTable(t_xsec_file, true);
}

inline void genie_set_event_generator_list(const std::string& t_evg) {
    genie::RunOpt::Instance()->SetEventGeneratorList(t_evg);
}

inline void genie_reconfiguration(bool ignore_alg_opt_out=false) {
    genie::AlgFactory::Instance()->ForceReconfiguration(ignore_alg_opt_out);
}

template<typename Value=double>
double E3_E0_p0_m3_m4_cosa(Value E0, Value p0, Value m3, Value m4, Value cosa)
{
    Value s = (E0 + p0) * (E0 - p0);
    Value sm32m42 = s + m3 * m3 - m4 * m4;
    Value p02cosa2 = p0 * p0 * cosa * cosa;
    Value E02p02cosa2 = E0 * E0 - p02cosa2;
    Value delta = p02cosa2 * sm32m42 * sm32m42 - 4.0 * E02p02cosa2 * p02cosa2 * m3 * m3;
    if (delta < 0) { return m3; }
    Value E3 = 0.0;
    if (cosa >= 0) {
        E3 = (E0 * sm32m42 + std::sqrt(delta)) / (2.0 * E02p02cosa2);
    } else {
        E3 = (E0 * sm32m42 - std::sqrt(delta)) / (2.0 * E02p02cosa2);
    }
    return std::max(E3, m3);
}


template<typename Value=double>
Value E3_E1_m1_m2_m3_m4_cosa(Value E1, Value m1, Value m2, Value m3, Value m4, Value cosa)
{
    Value E0 = E1 + m2;
    Value p0 = std::sqrt((E1 + m1) * (E1 - m1));
    return E3_E0_p0_m3_m4_cosa(E0, p0, m3, m4, cosa);
}

template<typename Value=double>
genie::Range1D_t Q2l_Echip(const genie::Interaction *interp, Value Echip,
                           Value Wcut=std::numeric_limits<Value>::max())
{
    const genie::InitialState & init_state = interp->InitState();
    double mN = init_state.Tgt().HitNucMass();
    double Echi = init_state.ProbeE(genie::kRfHitNucRest);
    auto Wl = interp->PhaseSpacePtr()->WLim();
    double Wmin = Wl.min;
    double Wmax = std::min(Wl.max, Wcut);
    double factor = mN*mN + 2*mN*(Echi - Echip);
    return genie::Range1D_t(factor - Wmax*Wmax, factor - Wmin*Wmin);
}

template<typename Value=double>
genie::Range1D_t EchipLim(const genie::Interaction *interp,
                          Value Wcut=std::numeric_limits<Value>::max(),
                          Value Q2cut=std::numeric_limits<Value>::max())
{
    const genie::InitialState & init_state = interp->InitState();
    double mN = init_state.Tgt().HitNucMass();
    double Echi = init_state.ProbeE(genie::kRfHitNucRest);
    auto Wl = interp->PhaseSpacePtr()->WLim();
    double Wmin = Wl.min;
    double Wmax = std::min(Wl.max, Wcut);
    interp->KinePtr()->SetW(Wmin);
    double Q2min = interp->PhaseSpacePtr()->Q2Lim_W().min;
    interp->KinePtr()->SetW(Wmax);
    double Q2max = interp->PhaseSpacePtr()->Q2Lim_W().max;
    Q2max = std::min(Q2max, Q2cut);

    double Echipmin = (2.0*Echi*mN + mN*mN - Q2max - Wmax*Wmax) / (2.0*mN);
    double Echipmax = (2.0*Echi*mN + mN*mN - Q2min - Wmin*Wmin) / (2.0*mN);
    double mchi = genie::PDGLibrary::Instance()->Find(genie::kPdgDarkMatter)->Mass();
    Echipmin = std::max(mchi, Echipmin);
    Echipmax = std::min(Echi, Echipmax);
    return genie::Range1D_t(Echipmin, Echipmax);
}


template<typename Value=double>
genie::Range1D_t EchipLim_ctl(const genie::Interaction *interp, Value ctl,
                              Value Wcut=std::numeric_limits<Value>::max(),
                              Value Q2cut=std::numeric_limits<Value>::max())
{
    auto Echipl = EchipLim(interp, Wcut, Q2cut);
    if (std::abs(ctl) >= 1) {
        return genie::Range1D_t(Echipl.max, Echipl.max);
    }
    double Echi = interp->InitState().ProbeE(genie::kRfHitNucRest);
    double mchi = interp->FSPrimLepton()->Mass();
    auto Wl = interp->PhaseSpacePtr()->WLim();
    double Wmin = Wl.min;
    double Wmax = std::min(Wl.max, Wcut);
    const genie::InitialState & init_state = interp->InitState();
    double mN = init_state.Tgt().HitNucMass();
    double Echipmin = E3_E1_m1_m2_m3_m4_cosa(Echi, mchi, mN, mchi, Wmax, ctl);
    double Echipmax = E3_E1_m1_m2_m3_m4_cosa(Echi, mchi, mN, mchi, Wmin, ctl);
    // std::cout << "Echipmin = " << Echipmin << " Echipmax = " << Echipmax << " Echipl.min" << Echipl.min << " Echipl.max = " << Echipl.max << " cos = " << ctl << std::endl;
    return genie::Range1D_t(std::max(Echipmin, Echipl.min),
                            std::min(Echipmax, Echipl.max));
}

template<typename Value=double>
genie::Range1D_t ctllimits_Echip(const genie::Interaction *interp, Value Echip,
                                 Value Wcut=std::numeric_limits<Value>::max(),
                                 Value Q2cut=std::numeric_limits<Value>::max())
{
    auto Echipl = EchipLim(interp, Wcut, Q2cut);
    if (Echipl.min >= Echipl.max || Echip <= Echipl.min || Echip >= Echipl.max) {
        return genie::Range1D_t(1.0, -1.0);
    }

    const genie::InitialState & init_state = interp->InitState();
    double mN = init_state.Tgt().HitNucMass();
    double Echi = init_state.ProbeE(genie::kRfHitNucRest);
    double mchi = interp->FSPrimLepton()->Mass();
    double Q2_W2 = mN * mN + 2.0 * mN * (Echi - Echip);
    auto Q2l = Q2l_Echip(interp, Echip, Wcut);
    Q2l.max = std::min(Q2l.max, Q2cut);
    if (Q2l.min >= Q2l.max) {
        return genie::Range1D_t(1.0, -1.0);
    }
    double denominator = std::sqrt((Echi+mchi)*(Echi-mchi)
            * (mN * (2.0*Echi - 2.0*mchi + mN) - Q2_W2)
            * (mN * (2.0*Echi + 2.0*mchi + mN) - Q2_W2));
    double numerator_factor = 2.0*Echi*Echi*mN - 2.0*mN*mchi*mchi
        + Echi*(mN*mN - Q2_W2);
    double ctlmin = (numerator_factor - mN * Q2l.max) / denominator;
    double ctlmax = (numerator_factor - mN * Q2l.min) / denominator;
    if (ctlmax <= -1) {
        return genie::Range1D_t(1.0, -1.0);
    }
    ctlmin = std::max(ctlmin, -1.0);
    ctlmax = std::min(ctlmax, 1.0);
    return genie::Range1D_t(ctlmin, ctlmax);
}

template<typename Value=double>
Value Q2_costheta(Value costheta, Value E1, Value m1, Value m2)
{
    Value ct2 = costheta * costheta;
    Value st2 = std::max(1.0 - ct2, 0.0);
    Value delta = ct2 * (m2*m2 - st2 * m1*m1);
    if (delta < 0) { return 0.0; }
    Value p02 = (E1 + m1) * (E1 - m1);
    Value denominator = (E1 + m2) * (E1 + m2) - p02 * ct2;
    Value Q2 = 0.0;
    if (costheta >= 0.0) {
        Q2 = 2.0 * m2 * p02 * (st2 * E1 + m2 - std::sqrt(delta)) / denominator;
    } else {
        Q2 = 2.0 * m2 * p02 * (st2 * E1 + m2 + std::sqrt(delta)) / denominator;
    }
    return std::max(Q2, 0.0);
}

template<typename Value=double>
Value dXSec_dQ2_DMCEL(Value Q2, Value Echi, int Z, int A)
{
    const genie::DMCELPXSec * xsec_cel
        = dynamic_cast<const genie::DMCELPXSec *>(
        genie::AlgFactory::Instance()->GetAlgorithm(
            "genie::DMCELPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMCEL(
            tgtc, genie::kPdgDarkMatter, Echi);

    interp->KinePtr()->SetQ2(Q2);
    double xsec = xsec_cel->XSec(interp, genie::kPSQ2fE);

    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dTchi_DMCEL(Value Tchi, Value Echi, int Z, int A)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    double mchi = pdglib->Find(genie::kPdgDarkMatter)->Mass();
    double mA = pdglib->Find(genie::pdg::IonPdgCode(A, Z))->Mass();
    double TA = Echi - Tchi - mchi;
    double Q2 = 2.0 * mA * TA;
    double xsec = dXSec_dQ2_DMCEL(Q2, Echi, Z, A);
    return xsec * 2.0 * mA;
}

template<typename Value=double>
Value dXSec_dQ2_DMEL_hitnuc(Value Q2, Value Echi, int Z, int A, int hitnuc)
{
    const genie::AhrensDMELPXSec * xsec_el
        = dynamic_cast<const genie::AhrensDMELPXSec *>(
                genie::AlgFactory::Instance()->GetAlgorithm(
                    "genie::AhrensDMELPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DME(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    interp->KinePtr()->SetQ2(Q2);
    double xsec = xsec_el->XSec(interp, genie::kPSQ2fE);

    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dQ2_DMEL(Value Q2, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dQ2_DMEL_hitnuc(Q2, Echi, Z, A, nucleon);
    }
    return xsec;
}


template<typename Value=double>
Value dXSec_dTchi_DMEL_hitnuc(Value Tchi, Value Echi, int Z, int A, int hitnuc)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    double mchi = pdglib->Find(genie::kPdgDarkMatter)->Mass();
    double mN = pdglib->Find(hitnuc)->Mass();
    double TN = Echi - Tchi - mchi;
    double Q2 = 2.0 * mN * TN;
    double xsec = dXSec_dQ2_DMEL_hitnuc(Q2, Echi, Z, A, hitnuc);
    return xsec * 2.0 * mN;
}

template<typename Value=double>
Value dXSec_dTchi_DMEL(Value Tchi, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dTchi_DMEL_hitnuc(Tchi, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dQ2_DMRES_hitnuc(Value Q2, Value Echi, int Z, int A, int hitnuc)
{
    constexpr std::array<genie::Resonance_t, 13> resonances = {
        genie::kP33_1232, genie::kS11_1535, genie::kD13_1520,
        genie::kS11_1650, genie::kD15_1675, genie::kS31_1620,
        genie::kD33_1700, genie::kP11_1440, genie::kP13_1720,
        genie::kF15_1680, genie::kP31_1910, genie::kF35_1905,
        genie::kF37_1950
    };
    const genie::DMRESPXSec* xsec_res
        = dynamic_cast<const genie::DMRESPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::DMRESPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMRES(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    double xsec = 0.0;
    for (auto the_res : resonances) {
        interp->ExclTagPtr()->SetResonance(the_res);
        auto Wl_res = interp->PhaseSpacePtr()->WLim();
        Wl_res.max = TMath::Min(Wl_res.max, xsec_res->Wcut(interp));
        genie::utils::gsl::dXSecInel_dQ2_E * func
            = new genie::utils::gsl::dXSecInel_dQ2_E(
                    xsec_res, interp, Wl_res.min, Wl_res.max);
        xsec += func->DoEval(Q2) * 1e-38 * genie::units::cm2;
        delete func;
    }
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dQ2_DMRES(Value Q2, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dQ2_DMRES_hitnuc(Q2, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dTchi_DMRES_hitnuc(Value Tchi, Value Echi, int Z, int A, int hitnuc)
{
    constexpr std::array<genie::Resonance_t, 13> resonances = {
        genie::kP33_1232, genie::kS11_1535, genie::kD13_1520,
        genie::kS11_1650, genie::kD15_1675, genie::kS31_1620,
        genie::kD33_1700, genie::kP11_1440, genie::kP13_1720,
        genie::kF15_1680, genie::kP31_1910, genie::kF35_1905,
        genie::kF37_1950
    };
    double mchi = genie::PDGLibrary::Instance()->Find(genie::kPdgDarkMatter)->Mass();
    const genie::DMRESPXSec* xsec_res = dynamic_cast<const genie::DMRESPXSec *>(
        genie::AlgFactory::Instance()->GetAlgorithm("genie::DMRESPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMRES(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    double xsec = 0.0;
    for (auto the_res : resonances) {
        interp->ExclTagPtr()->SetResonance(the_res);
        auto ctll = ctllimits_Echip(
                interp, Tchi+mchi, xsec_res->Wcut(interp), xsec_res->Q2cut());
        if (ctll.min < ctll.max) {
            genie::utils::gsl::dXSecInel_dTl_E * func =
                new genie::utils::gsl::dXSecInel_dTl_E(
                    xsec_res, interp, ctll.min, ctll.max);
            xsec += func->DoEval(Tchi) * 1e-38 * genie::units::cm2;
            delete func;
        }
    }
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dTchi_DMRES(Value Tchi, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dTchi_DMRES_hitnuc(Tchi, Echi, Z, A, nucleon);
    }
    return xsec;
}


template<typename Value=double>
Value dXSec_dQ2_DMDIS_hitnuc(Value Q2, Value Echi, int Z, int A, int hitnuc)
{
    const genie::QPMDMDISPXSec* xsec_dis
        = dynamic_cast<const genie::QPMDMDISPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::QPMDMDISPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMDI(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    auto Wl_dis = interp->PhaseSpacePtr()->WLim();
    genie::utils::gsl::dXSecInel_dQ2_E * func =
        new genie::utils::gsl::dXSecInel_dQ2_E(
                xsec_dis, interp, Wl_dis.min, Wl_dis.max);
    double xsec = func->DoEval(Q2) * 1e-38 * genie::units::cm2;
    delete func;
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dQ2_DMDIS(Value Q2, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dQ2_DMDIS_hitnuc(Q2, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dTchi_DMDIS_hitnuc(Value Tchi, Value Echi, int Z, int A, int hitnuc)
{
    const genie::QPMDMDISPXSec* xsec_dis
        = dynamic_cast<const genie::QPMDMDISPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::QPMDMDISPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMDI(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    double mchi = genie::PDGLibrary::Instance()->Find(genie::kPdgDarkMatter)->Mass();
    auto ctll = ctllimits_Echip(interp, Tchi+mchi);
    genie::utils::gsl::dXSecInel_dTl_E * func = new genie::utils::gsl::dXSecInel_dTl_E(
                xsec_dis, interp, ctll.min, ctll.max);
    double xsec = func->DoEval(Tchi) * 1e-38 * genie::units::cm2;
    delete func;
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dTchi_DMDIS(Value Q2, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dTchi_DMDIS_hitnuc(Q2, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dcostheta_DMCEL(Value ctl, Value Echi, int Z, int A)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    double mchi = pdglib->Find(genie::kPdgDarkMatter)->Mass();
    double mA = pdglib->Find(genie::pdg::IonPdgCode(A, Z))->Mass();
    double Q2 = Q2_costheta(ctl, Echi, mchi, mA);
    double k = std::sqrt((Echi + mchi) * (Echi - mchi));
    double Echip = Echi - Q2 / (2.0 * mA);
    double kp = std::sqrt((Echip + mchi) * (Echip - mchi));
    double xsec = dXSec_dQ2_DMCEL(Q2, Echi, Z, A);
    return xsec * 2.0 * k * kp;
}

template<typename Value=double>
Value dXSec_dcostheta_DMEL_hitnuc(Value ctl, Value Echi, int Z, int A, int hitnuc)
{
    genie::PDGLibrary * pdglib = genie::PDGLibrary::Instance();
    double mchi = pdglib->Find(genie::kPdgDarkMatter)->Mass();
    double mN = pdglib->Find(hitnuc)->Mass();
    double Q2 = Q2_costheta(ctl, Echi, mchi, mN);
    double k = std::sqrt((Echi + mchi) * (Echi - mchi));
    double Echip = Echi - Q2 / (2.0 * mN);
    double kp = std::sqrt((Echip + mchi) * (Echip - mchi));
    double xsec = dXSec_dQ2_DMEL_hitnuc(Q2, Echi, Z, A, hitnuc);
    return xsec * 2.0 * k * kp;
}

template<typename Value=double>
Value dXSec_dcostheta_DMEL(Value ctl, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dcostheta_DMEL_hitnuc(ctl, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dcostheta_DMRES_hitnuc(Value ctl, Value Echi, int Z, int A, int hitnuc)
{
    constexpr std::array<genie::Resonance_t, 13> resonances = {
        genie::kP33_1232, genie::kS11_1535, genie::kD13_1520,
        genie::kS11_1650, genie::kD15_1675, genie::kS31_1620,
        genie::kD33_1700, genie::kP11_1440, genie::kP13_1720,
        genie::kF15_1680, genie::kP31_1910, genie::kF35_1905,
        genie::kF37_1950
    };
    const genie::DMRESPXSec* xsec_res = dynamic_cast<const genie::DMRESPXSec *>(
        genie::AlgFactory::Instance()->GetAlgorithm("genie::DMRESPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMRES(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    double xsec = 0.0;
    for (auto the_res : resonances) {
        interp->ExclTagPtr()->SetResonance(the_res);
        auto Tll = EchipLim_ctl(interp, ctl, xsec_res->Wcut(interp), xsec_res->Q2cut());
        double mchi = genie::PDGLibrary::Instance()->Find(genie::kPdgDarkMatter)->Mass();
        Tll.min -= mchi;
        Tll.max -= mchi;
        if (Tll.min < Tll.max) {
            genie::utils::gsl::dXSecInel_dctl_E * func =
                new genie::utils::gsl::dXSecInel_dctl_E(
                    xsec_res, interp, Tll.min, Tll.max);
            xsec += func->DoEval(ctl) * 1e-38 * genie::units::cm2;
            delete func;
        }
    }
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dcostheta_DMRES(Value ctl, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dcostheta_DMRES_hitnuc(ctl, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value dXSec_dcostheta_DMDIS_hitnuc(Value ctl, Value Echi, int Z, int A, int hitnuc)
{
    const genie::QPMDMDISPXSec* xsec_dis
        = dynamic_cast<const genie::QPMDMDISPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::QPMDMDISPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMDI(
        tgtc, hitnuc, genie::kPdgDarkMatter, Echi);

    double mchi = genie::PDGLibrary::Instance()->Find(genie::kPdgDarkMatter)->Mass();
    auto Tll = EchipLim_ctl(interp, ctl);
    Tll.min -= mchi;
    Tll.max -= mchi;
    genie::utils::gsl::dXSecInel_dctl_E * func = new genie::utils::gsl::dXSecInel_dctl_E(
                xsec_dis, interp, Tll.min, Tll.max);
    double xsec = func->DoEval(ctl) * 1e-38 * genie::units::cm2;
    delete func;
    delete interp;
    return xsec;
}

template<typename Value=double>
Value dXSec_dcostheta_DMDIS(Value ctl, Value Echi, int Z, int A)
{
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        xsec += dXSec_dcostheta_DMDIS_hitnuc(ctl, Echi, Z, A, nucleon);
    }
    return xsec;
}

template<typename Value=double>
Value total_xsec_DMCEL(Value Echi, int Z, int A)
{
    const genie::DMCELPXSec * xsec_cel
        = dynamic_cast<const genie::DMCELPXSec *>(
        genie::AlgFactory::Instance()->GetAlgorithm(
            "genie::DMCELPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMCEL(
            tgtc, genie::kPdgDarkMatter, Echi);

    double xsec = xsec_cel->Integral(interp);

    delete interp;
    return xsec;
}

template<typename Value=double>
Value total_xsec_DMEL(Value Echi, int Z, int A)
{
    const genie::AhrensDMELPXSec * xsec_el
        = dynamic_cast<const genie::AhrensDMELPXSec *>(
                genie::AlgFactory::Instance()->GetAlgorithm(
                    "genie::AhrensDMELPXSec", "Velocity0"));


    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (int nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        const genie::Interaction * interp = genie::Interaction::DME(
            tgtc, nucleon, genie::kPdgDarkMatter, Echi);
        xsec += xsec_el->Integral(interp);
        delete interp;
    }

    return xsec;
}

template<typename Value=double>
Value total_xsec_DMRES(Value Echi, int Z, int A)
{
    constexpr std::array<genie::Resonance_t, 13> resonances = {
        genie::kP33_1232, genie::kS11_1535, genie::kD13_1520,
        genie::kS11_1650, genie::kD15_1675, genie::kS31_1620,
        genie::kD33_1700, genie::kP11_1440, genie::kP13_1720,
        genie::kF15_1680, genie::kP31_1910, genie::kF35_1905,
        genie::kF37_1950
    };
    const genie::DMRESPXSec* xsec_res
        = dynamic_cast<const genie::DMRESPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::DMRESPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);

    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (auto nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }
        const genie::Interaction * interp = genie::Interaction::DMRES(
            tgtc, nucleon, genie::kPdgDarkMatter, Echi);
        for (auto the_res : resonances) {
            interp->ExclTagPtr()->SetResonance(the_res);
            xsec += xsec_res->Integral(interp);
        }
        delete interp;
    }

    return xsec;
}

template<typename Value=double>
Value total_xsec_DMDIS(Value Echi, int Z, int A)
{
    const genie::QPMDMDISPXSec* xsec_dis
        = dynamic_cast<const genie::QPMDMDISPXSec *>(
            genie::AlgFactory::Instance()->GetAlgorithm(
                "genie::QPMDMDISPXSec", "Velocity0"));

    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const std::array<int, 2> nucleons {genie::kPdgProton, genie::kPdgNeutron};
    double xsec = 0;
    for (auto nucleon : nucleons) {
        if (nucleon == genie::kPdgProton && Z <= 0) { continue; }
        if (nucleon == genie::kPdgNeutron && A - Z <= 0) { continue; }

        const genie::Interaction * interp = genie::Interaction::DMDI(
            tgtc, nucleon, genie::kPdgDarkMatter, Echi);
        xsec += xsec_dis->Integral(interp);
        delete interp;
    }
    return xsec;
}

template<typename Value=double>
std::vector<Value> limits_Echip_DMRES(Value Echi, int Z, int A,
                        Value Wcut, Value Q2cut)
{
    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMRES(
        tgtc, genie::kPdgProton, genie::kPdgDarkMatter, Echi);

    auto Elim = EchipLim(interp, Wcut, Q2cut);
    delete interp;
    return {Elim.min, Elim.max};
}

template<typename Value=double>
std::vector<double> limits_Echip_DMDIS(Value Echi, int Z, int A)
{
    int tgtc = genie::pdg::IonPdgCode(A, Z);
    const genie::Interaction * interp = genie::Interaction::DMDI(
        tgtc, genie::kPdgProton, genie::kPdgDarkMatter, Echi);

    auto Elim = EchipLim<Value>(interp);
    delete interp;
    return {Elim.min, Elim.max};
}

template<typename Value=double>
std::vector<Value> limits_Q2_W(Value Echi, Value M, Value mchi, Value W,
                               Value Q2min_cut=genie::controls::kMinQ2Limit_DM)
{
    auto Q2lim = genie::utils::kinematics::DarkQ2Lim_W(Echi, M, mchi, W, Q2min_cut);
    return {Q2lim.min, Q2lim.max};
}

template<typename Value=double>
std::vector<Value> limits_Q2(Value Echi, Value M, Value mchi,
                             Value Q2min_cut=genie::controls::kMinQ2Limit_DM)
{
    auto Q2lim = genie::utils::kinematics::DarkQ2Lim(Echi, M, mchi, Q2min_cut);
    return {Q2lim.min, Q2lim.max};
}

template<typename Value=double>
std::vector<Value> limits_W(Value Echi, Value M, Value mchi)
{
    auto Wlim = genie::utils::kinematics::DarkWLim(Echi, M, mchi);
    return {Wlim.min, Wlim.max};
}

template<typename Value=double>
std::vector<Value> limits_costheta_Echip(Value Echi, Value Echip)
{
    int tgtc = genie::pdg::IonPdgCode(1, 1);
    const genie::Interaction * interp = genie::Interaction::DMDI(
        tgtc, genie::kPdgProton, genie::kPdgDarkMatter, Echi);
    auto ctllim = ctllimits_Echip(interp, Echip);
    delete interp;
    return {ctllim.min, ctllim.max};
}

} // namespace darkprop
