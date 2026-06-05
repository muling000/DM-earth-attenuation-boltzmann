/**
 * @file
 * Wrapper for the boosted DM model implemented in the GENIE code.
 */
#pragma once
#include <map>
#include <cmath>
#include <TBits.h>
#include <TLorentzVector.h>
#include <Framework/Algorithm/AlgConfigPool.h>
#include <Framework/EventGen/GEVGDriver.h>
#include <Framework/EventGen/EventRecord.h>
#include <Framework/Interaction/Interaction.h>
#include <Framework/ParticleData/PDGCodes.h>
#include <Framework/ParticleData/PDGLibrary.h>
#include <Framework/ParticleData/PDGUtils.h>
#include <Framework/GHEP/GHepParticle.h>
#include <Framework/Utils/RunOpt.h>
#include <Framework/Utils/AppInit.h>
#include <Framework/Utils/XSecSplineList.h>
#include <Framework/Numerical/Spline.h>
#include <spdlog/spdlog.h>
#include "Base.hpp"
#include "Utils.hpp"
#include "Vector3d.hpp"
#include "GENIEWrapper.hpp"


namespace darkprop {

template<typename Vector3=Vector3d<double>, typename Value=double>
class GENIEBDM : public Particle<Vector3, Value>
{
    using type = GENIEBDM<Vector3, Value>;
private:
    // We use a simple map instead of the Pool provided by GENIE.
    std::map<int, genie::GEVGDriver*> drivers;
public:
    std::string event_gen_list;
    Value Tmin;
    Value Tmax;
    Value mmed;
    Value gz_org;
    Value gz;
public:
    explicit GENIEBDM(Value t_sigma0=1e-30*units::cm2,
                      Value t_t=0.0,
                      Value t_Tmin=1.1*units::GeV,
                      Value t_Tmax=10.0*units::GeV,
                      Value t_gz_org=1.0) {
        genie::PDGLibrary * const pdglib = genie::PDGLibrary::Instance();
        const TParticlePDG * const dmp = pdglib->Find(genie::kPdgDarkMatter);
        if (!dmp) { throw std::runtime_error("Dark matter is not added into GENIE."); }
        this->m = dmp->Mass();
        this->mmed = pdglib->Find(genie::kPdgMediator)->Mass();
        this->sigma0 = t_sigma0;
        this->t = t_t;
        this->Tmin = t_Tmin;
        this->Tmax = t_Tmax;
        this->gz_org = t_gz_org;
        this->in_medium = false;
        setgz();
    }
    ~GENIEBDM() { for (auto& [key, driver] : this->drivers) { delete driver; } }

    Value totalCrossSection(const Target& t) const override {
        int pdgt = genie::pdg::IonPdgCode(t.A, t.Z);
        Value xsec = const_cast<genie::Spline*>(
                const_cast<type*>(this)->drivers[pdgt]->XSecSumSpline()
            )->Evaluate(this->T + this->m);
        return xsec * std::pow(this->gz / this->gz_org, 4.0);
    }
    Value totalCrossSection(const std::string& splkey);

    Value scatter(const Target& t, RandomNumber<Value>& rn) override;

    void initGEVGDriver(int Z, int A, std::size_t knots=1000);
    void initGEVGDriver(const Medium<Vector3, Value>& m, std::size_t knots=1000) {
        spdlog::info("initiating GEVGDriver");
        for (const auto& t : m.getTargets()) {
            initGEVGDriver(t.Z, t.A, knots);
        }
    }

    void setgz(Value t_gz) { this->gz = t_gz; }
    void setgz() {
        this->gz = this->mmed * std::pow(M_PI * this->sigma0
            / (9.0 * std::pow(reduce_m(this->m, constants::mp), 2.0)), 1.0/4.0);
    }
};

template<typename Vector3, typename Value>
void GENIEBDM<Vector3, Value>::initGEVGDriver(int Z, int A, std::size_t knots)
{
    spdlog::info("init Z = {0:d}, A = {1:d}", Z, A);
    int pdgt = genie::pdg::IonPdgCode(A, Z);
    if (this->drivers.count(pdgt) != 0) { delete drivers[pdgt]; }
    this->drivers[pdgt] = new genie::GEVGDriver();
    this->drivers[pdgt]->SetEventGeneratorList(
        genie::RunOpt::Instance()->EventGeneratorList());
    this->drivers[pdgt]->Configure(genie::kPdgDarkMatter, Z, A);
    this->drivers[pdgt]->CreateSplines(knots, 1.1 * this->m + this->Tmax, true);
    // Although CreateXSecSumSpline assuming massless probe, we always
    // ensure that the splines have been created by CreateSplines so that
    // only nup4.Energy() is called by those splines, which is intact with
    // a massive DM probe.
    this->drivers[pdgt]->CreateXSecSumSpline(
        knots, this->m / 2, 1.1 * this->m + this->Tmax, true);
    this->drivers[pdgt]->SetUnphysEventMask(*genie::RunOpt::Instance()->UnphysEventMask());
}

template<typename Vector3, typename Value>
Value GENIEBDM<Vector3, Value>::scatter(const Target& t, RandomNumber<Value>& rn)
{
    int pdgt = genie::pdg::IonPdgCode(t.A, t.Z);
    const auto & p3 = this->p3;
    // px, py, pz, E (GeV)
    TLorentzVector dm_p4(p3.x, p3.y, p3.z, this->T + this->m);
    genie::EventRecord * event = this->drivers[pdgt]->GenerateEvent(dm_p4);
    if (!event) {
        dm_p4.Print();
        spdlog::error("EvGenList: {0}, gz = {1:.3f}, mchi = {2:.3e}, mmed = {3:.3e}",
            genie::RunOpt::Instance()->EventGeneratorList(), this->gz, this->m, this->mmed);
        throw std::runtime_error("!! GenerateEvent failed!");
    }

    const genie::GHepParticle *dm = event->FinalStatePrimaryLepton();
    // get DM final 4-momenta (@ LAB)
    if (dm) {
        const TLorentzVector& p4 = *(dm->P4());
        this->setP(Vector3(p4(0), p4(1), p4(2)));
    } else {
        spdlog::error("Final state DM does not exist!");
        throw std::runtime_error("Final state DM does not exist!");
    }

    delete event;
    return 1.0;
}

template<typename Vector3, typename Value>
Value GENIEBDM<Vector3, Value>::totalCrossSection(const std::string& splkey)
{
    genie::XSecSplineList * xssl = genie::XSecSplineList::Instance();
    bool spline_exists = xssl->SplineExists(splkey);
    if (!spline_exists) {
        spdlog::warn("{} doesn't exist! 0 returned.", splkey);
        return 0.0;
    }
    return xssl->GetSpline(splkey)->Evaluate(this->T + this->m);
}

} // namespace darkprop
