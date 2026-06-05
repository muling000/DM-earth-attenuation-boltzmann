/**
 * @file
 * Solar DM model. Scalar DM couples to proton.
 */

#pragma once

#include "Base.hpp"
#include "Sun.hpp"
#include "Vector3d.hpp"
#include "TwoBodyPhaseSpace.hpp"

namespace darkprop {

/**
 * The Solar DM model.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SolarDM : public Particle<Vector3, Value>
{
private:
    Value temp_scale = 1.0;
public:
    explicit SolarDM(bool t_in_medium=true) : Particle<Vector3, Value> {t_in_medium} {};
    explicit SolarDM(Value t_sigma0=1e-30*units::cm2,
                     Value t_m=1.*units::keV,
                     Value t_t=0.0) {
                     this->t = t_t;
                     this->sigma0 = t_sigma0;
                     this->m = t_m;
                     this->in_medium = true;
    }
    virtual ~SolarDM() {}

    /**
     * Set the temperature scaling factor multiplied to the real temperature.
     *
     * This is used to adjust the temperature of the Sun.
     */
    void setTempScale(Value t_temp_scale) {
        if (t_temp_scale <= 0) { throw std::runtime_error("temperature <= 0"); }
        this->temp_scale = t_temp_scale;
    }
    Value getTempScale() { return this->temp_scale; }

    Value totalCrossSection(const Target& t) const override {
        return t.Z * t.Z * this->sigma0;
    }
    Value scatter(const Target& t, RandomNumber<Value>& rn) override;
};

template<typename Vector3, typename Value>
Value SolarDM<Vector3, Value>::scatter(const Target& t, RandomNumber<Value>& rn)
{
    Value temperature = solar_temperature(this->r.norm())*this->temp_scale;
    Value s = std::sqrt(t.m * temperature);  // the standard variance sigma
    Vector3 pt(rn.normal(0.0, s), rn.normal(0.0, s), rn.normal(0.0, s));
    Value Et = std::sqrt(pt.dot(pt) + t.m * t.m);

    const Vector3& p1in = this->p3;

    Value E0 = Et + this->T + this->m;
    Vector3 p0 = pt + p1in;
    Value p0_norm = p0.norm();

    const auto [E1, phi1] = twobodyphase::sample_uniform(E0, p0_norm, this->m, t.m, rn);
    Value costheta1 = twobodyphase::costheta1(E0, p0_norm, this->m, t.m, E1);
    Value p1out = std::sqrt(E1 * E1 - this->m * this->m);

    this->setP(twobodyphase::final_momentum(p0, p1in, p1out, costheta1, phi1));
    return 1.0;
}

} // namespace darkprop
