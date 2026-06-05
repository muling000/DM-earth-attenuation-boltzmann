/**
 * @file
 * Base class for elastic scattering particle model.
 */

#pragma once

#include <spdlog/spdlog.h>
#include "Base.hpp"
#include "Utils.hpp"

namespace darkprop {


/**
 * Base class for elastic scattering particle model.
 */
template<typename Vector3, typename Value>
class ParticleElastic : public Particle<Vector3, Value>
{
    using base = Particle<Vector3, Value>;
public:
    ParticleElastic(Value t_sigma0=1e-30*units::cm*units::cm,
         Value t_m=1.*units::GeV,
         Value t_t=0.0) {
        this->t = t_t;
        this->sigma0 = t_sigma0;
        this->m = t_m;
        this->in_medium = false;
    }
    ParticleElastic(bool t_in_medium=false) : base {t_in_medium} {};
    virtual ~ParticleElastic() {}

    /**
     * Sample recoil kinetic energy of the target particle in laboratory frame.
     */
    virtual Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) = 0;

    /**
     * Return cosine of the scattering angle of DM particle in the target rest frame.
     * @param t the target
     * @param Tr recoil kinetic energy of the target
     */
    virtual Value scatteringCosTheta(const Target& t, Value Tr);

    /**
     * Sample final state DM particle using ParticleElastic::scatteringSampleTr and
     * ParticleElastic::scatteringCosTheta.
     */
    virtual Value scatter(const Target& t, RandomNumber<Value>& rn) override {
        Value Tr = scatteringSampleTr(t, rn);
        Value costheta = scatteringCosTheta(t, Tr);
        scatter_fixed_target(*this, this->T - Tr, costheta, rn);
        return 1.0;
    }

    /**
     * The maximal recoil kinetic energy of a fixed \p target in elastic scattering.
     */
    virtual Value maximumRecoilT(const Target& target) const {
        return this->T * (this->T + 2. * this->m)
               / (this->T + std::pow(this->m + target.m, 2) / (2. * target.m));
    }
};

template<typename Vector3, typename Value>
Value ParticleElastic<Vector3, Value>::scatteringCosTheta(const Target& t, Value Tr)
{
    // it seems no way to ensure |costheta| <= 1 :(
    Value p12 = this->T * (this->T + 2. * this->m);  // p1 squared
    Value Tf = this->T - Tr;  // final energy of DM
    // Value costheta = (p12 - Tr * (this->T + this->m + t.m))
    //                  / std::sqrt(p12 * Tf * (Tf + 2. * this->m));

    Value numerator = p12 - Tr * (this->T + this->m + t.m);
    Value costheta = std::sqrt(std::pow(numerator, 2.0)
                        / p12 / (Tf * (Tf + 2. * this->m)));
    costheta = costheta * (numerator > 0 ? 1 : -1);

    if (costheta > 1.0) {
        if (costheta > 1.0 + 1e-15) {
            spdlog::warn("costheta = {0:.18e}, where mchi = {1:.18e}, T = {2:.18e} "
                         "Tr = {3:.18e}, Tf = {4:.18e}. set costheta = 1",
                         costheta, this->m, this->T, Tr, Tf);
        }
        costheta = 1.0;
    } else if (costheta < -1.0) {
        if (costheta < -(1.0 + 1e-15)) {
            spdlog::warn("costheta = {0:.18e}, where mchi = {1:.18e}, T = {2:.18e} "
                         "Tr = {3:.18e}, Tf = {4:.18e}. set costheta = -1",
                         costheta, this->m, this->T, Tr, Tf);
        }
        costheta = -1.0;
    }
    return costheta;
}

} // namespace darkprop
