/**
 * @file
 * DM-electron interaction model. Currently only vector mediator is implemented.
 */

#pragma once

#include <stdexcept>
#include <spdlog/spdlog.h>
#include "Const.hpp"
#include "Vector3d.hpp"
#include "ParticleElastic.hpp"
#include "numerical/RootFinding.hpp"

namespace darkprop {

/**
 * A DM-electron vector coupling model. Used to cross check for arXiv:2403.08361.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class DMElectron : public ParticleElastic<Vector3, Value>
{
private:
    using Base = ParticleElastic<Vector3, Value>;
    static constexpr Value me = constants::me;
    static constexpr Value alpha2 = constants::alpha2;
    Value mr;
    Value mmed;
    Value A() const { return 2.0 * me * std::pow(this->m + this->T, 2.0); }
    Value B() const { return -(2.0*me*this->T + std::pow(this->m+me, 2.0)); }
    Value C() const { return me; }
public:
    std::string mediator_limit = "";
    Value abstol = 0;
    Value reltol = 1e-12;
    std::size_t maxit = 100;
    explicit DMElectron(Value t_sigma0=1e-30*units::cm*units::cm,
                        Value t_m=1.*units::GeV,
                        Value t_mmed=1.*units::GeV,
                        Value t_t=0.0) : Base(t_sigma0, t_m, t_t) {
        this->mmed = t_mmed;
        this->mr = this->mmed * this->mmed / (2.0 * me);
    }

    Value FDM(Value Te) const {
        Value poly = A() + B() * (Te + C() * Te) / (2.0*me*this->m*this->m);
        if (mediator_limit == "heavy") {
            return poly;
        }
        if (this->mmed <= 0) {
            throw std::invalid_argument("mediator mass <= 0");
        }
        return std::pow(mr+alpha2*me/2.0, 2.0) / std::pow(mr + Te, 2.0) * poly;
    }
    Value FDMintegral(Value Tp) const {
        Value a = A();
        Value b = B();
        Value c = C();
        Value fac = (2.0 * me * this->m * this->m);
        if (mediator_limit == "heavy") {
            Value intg = a * (Tp + b / 2.0 * (Tp + c / 3.0 * Tp)) / fac;
            if (!(intg >= 0)) {
                spdlog::error("FDMintegral: {0:.16e} T = {1:.16e}", intg, this->T);
                throw std::runtime_error("FDMintegral error");
            }
            return intg;
        }
        if (this->mmed <= 0) {
            throw std::invalid_argument("mediator mass <= 0");
        }
        Value intg = std::pow(mr + alpha2 * me / 2.0, 2.0) / fac * (
            (a / mr - b + c * (Tp + 2.0 * mr)) * Tp / (Tp + mr)
             + (2.0 * mr * c - b) * std::log(mr / (Tp + mr)));
        // Value x = mr / (Tp + mr);
        // Value logx = std::log(x);
        // Value intg = std::pow(mr + alpha2 * me / 2.0, 2.0) / fac * (
        //         a * x / mr - b * (x + logx)
        //         + c * ((Tp + 2.0*mr) * x + 2.0 * mr * logx));
        if (!(intg > 0)) {
            spdlog::error("FDMintegral: {0:.16e} T = {1:.16e}", intg, this->T);
            throw std::runtime_error("FDMintegral error");
        }
        return intg;
    }

    Value inverseCDF(Value xi, const Target& t) {
        Value Trmax = this->maximumRecoilT(t);
        Value xi_fdm_integral_max = xi * FDMintegral(Trmax);
        return numerical::rtsafe<Value>(
            [=](Value Te) {
                return this->FDMintegral(Te) - xi_fdm_integral_max;
            },
            [=](Value Te) {
                return this->FDM(Te);
            },
            0.0, Trmax, this->abstol, this->reltol, this->maxit);
    }

    Value totalCrossSection(const Target& t) const override {
        if (t.A != -1) {
            throw(std::invalid_argument("target is not electron"));
        }
        Value Trmax = this->maximumRecoilT(t);
        return this->sigma0 / Trmax * std::pow(this->m + me, 2.0)
            / (2.0 * me * this->T + std::pow(this->m + me, 2.0))
            * FDMintegral(Trmax);
    }

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override
    {
        Value xi = rn.uniform_xi();
        return inverseCDF(xi, t);
    }
};

} // namespace darkprop
