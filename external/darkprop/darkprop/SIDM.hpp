/**
 * @file
 * Spin independent elastic scattering DM model with constant cross section.
 */

#pragma once

#include "ParticleElastic.hpp"
#include "Const.hpp"
#include "Vector3d.hpp"

namespace darkprop {

/**
 * Spin independent elastic scattering DM model with constant cross section.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SIDM : public ParticleElastic<Vector3, Value>
{
    using Base = ParticleElastic<Vector3, Value>;
public:
    // TODO Base class constructor
    explicit SIDM(Value t_sigma0=1e-30*units::cm*units::cm,
                  Value t_m=1.*units::GeV,
                  Value t_t=0.0) : Base(t_sigma0, t_m, t_t) { }
    SIDM(Value t_sigma0,
         Value t_m,
         Value t_t,
         Vector3 t_r,
         Vector3 t_v)
    {
        this->sigma0 = t_sigma0;
        this->m = t_m;
        this->t = t_t;
        this->r = t_r;
        this->v = t_v;
        this->updateP3TwithV();
        this->updateEp();
        this->updateEr();
        this->in_medium = false;
    }
    SIDM(Vector3 t_r,
         Vector3 t_p,
         Value t_sigma0,
         Value t_m,
         Value t_t)
    {
        this->sigma0 = t_sigma0;
        this->m = t_m;
        this->t = t_t;
        this->r = t_r;
        this->p3 = t_p;
        this->T = this->TFromP3();
        this->updateVwithP3T();
        this->updateEp();
        this->updateEr();
        this->in_medium = false;
    }
    ~SIDM() {}
    Value sigmaEff(const Target& t) const {
        if (t.A == -1) {
            // for electron
            return this->sigma0;
        } else {
            return this->sigma0 * t.A*t.A * std::pow(t.m / constants::mu *
                        (this->m + constants::mu) / (this->m + t.m), 2.0);
        }
    }
    Value differentialCrossSection(Value Tr, const Target& t) const {
        Value Trmax = this->maximumRecoilT(t);
        if (Tr > Trmax || Tr < 0) return 0.0;
        return sigmaEff(t) / Trmax;
    }
    Value totalCrossSection(const Target& t) const override {
        return sigmaEff(t);
    }
    // Value sampleTanTheta(const Target& t, RandomNumber<Value>& rn) {
    //     Value xi = rn.uniform_xi();

    //     Value coschi = 1. - 2.*xi;
    //     Value sinchi = std::sqrt(1. - coschi*coschi);
    //     Value m1plusm2 = this->m + t.m;
    //     Value E1plusm2 = this->T + m1plusm2;

    //     return sinchi
    //         / (E1plusm2/std::sqrt(2.*t.m*this->T+m1plusm2*m1plusm2)) //gamma
    //         / (coschi + (this->m*m1plusm2+this->T*t.m)/(t.m * E1plusm2));
    // }

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override
    {
        return this->maximumRecoilT(t) * rn.uniform_xi();
    }
};

} // namespace darkprop
