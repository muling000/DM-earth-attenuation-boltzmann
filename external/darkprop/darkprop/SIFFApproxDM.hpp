#include "SIDM.hpp"
#include "Const.hpp"

namespace darkprop {

/**
 * SIDM with exponential form factor.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SIFFExpDM : public SIDM<Vector3, Value>
{
private:
    const Value s2 = 1.0 * units::fm * units::fm;
public:
    SIFFExpDM() : SIDM<Vector3, Value>() {}
    ~SIFFExpDM() {}

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override
    {
        Value q2max = 2.0 * t.m * this->maximumRecoilT(t);
        Value xi = rn.uniform_xi();
        Value q2 = -std::log(1.0-(1.0-std::exp(-q2max*s2))*xi) / s2;
        return q2 / (2.0*t.m);
    }

    Value totalCrossSection(const Target& t) const override {
        Value q2max = 2.0 * t.m * this->maximumRecoilT(t);
        return this->sigmaEff(t) * (1-std::exp(-q2max*s2)) / s2 / q2max;
    }
};

/**
 * SIDM with Heaviside step function form factor.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SIFFCutDM : public SIDM<Vector3, Value>
{
public:
    Value q2cutoff;
    SIFFCutDM(Value t_q2cutoff=2.0/(units::fm*units::fm))
        : SIDM<Vector3, Value>(), q2cutoff {t_q2cutoff} {}
    ~SIFFCutDM() {}

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override
    {
        Value q2max = 2.0 * t.m * this->maximumRecoilT(t);
        Value q2 = std::min(q2cutoff, q2max) * rn.uniform_xi();
        return q2 / (2.0*t.m);
    }

    Value totalCrossSection(const Target& t) const override {
        Value q2max = 2.0 * t.m * this->maximumRecoilT(t);
        Value ff2int = std::min(q2cutoff, q2max);
        return this->sigmaEff(t) * ff2int / q2max;
    }
};

} // namespace darkprop
