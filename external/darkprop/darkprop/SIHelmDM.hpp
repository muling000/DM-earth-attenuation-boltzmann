/**
 * @file
 * Spin independent elastic scattering DM model with Helm form factor.
 */

#pragma once

#include <map>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "SIDM.hpp"
#include "Utils.hpp"
#include "Vector3d.hpp"
#include "FormFactor.hpp"

namespace darkprop {

/**
 * Spin independent elastic scattering DM model with Helm form factor.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SIHelmDM : public SIDM<Vector3, Value>
{
private:
    gsl_integration_workspace *w;
    Value Tr_min_init = 0;  // store lower bound
    Value Tr_max_init = 0;  // store upper bound
    std::map<std::string, gsl_spline*> map_cdf;
    std::map<std::string, gsl_spline*> map_cdf_inv;
    std::map<std::string, gsl_interp_accel*> map_acc;
    std::map<std::string, gsl_interp_accel*> map_acc_inv;
    void alloc_gsl() {
        w = gsl_integration_workspace_alloc(1000);
    }
public:
    explicit SIHelmDM(Value t_sigma0=1e-30*units::cm*units::cm,
                      Value t_m=1.*units::GeV,
                      Value t_t=0.0)
        : SIDM<Vector3, Value>(t_sigma0, t_m, t_t) {
        alloc_gsl();
        this->in_medium = false;
    }

    SIHelmDM(Value t_sigma0, Value t_m, Value t_t, Vector3 t_r, Vector3 t_v)
        : SIDM<Vector3, Value>(t_sigma0, t_m, t_t, t_r, t_v) {
        alloc_gsl();
        this->in_medium = false;
    }

    ~SIHelmDM() {
        gsl_integration_workspace_free(w);
        for (auto it = map_cdf.begin(); it != map_cdf.end(); ++it) {
            gsl_spline_free(it->second);
        }
        for (auto it = map_cdf_inv.begin(); it != map_cdf_inv.end(); ++it) {
            gsl_spline_free(it->second);
        }
        for (auto it = map_acc.begin(); it != map_acc.end(); ++it) {
            gsl_interp_accel_free(it->second);
        }
        for (auto it = map_acc_inv.begin(); it != map_acc_inv.end(); ++it) {
            gsl_interp_accel_free(it->second);
        }
    }

    SIHelmDM& operator=(const SIHelmDM& other) = delete;

    Value differentialCrossSection(Value Tr, const Target& t) const {
        Value ff;
        if (t.A >= 7) {
            ff = ff_Helm(2. * t.m * Tr, t.A);
        } else {
            ff = ff_dipole(2. * t.m * Tr, t.A);
        }
        return SIDM<Vector3, Value>::differentialCrossSection(Tr, t)
            * std::pow(ff, 2.0);
    }

    static double integrand(double q2, void *t_A) {
        int *A = static_cast<int *>(t_A);
        if (*A >= 7) {
            return std::pow(ff_Helm(q2, *A), 2.);
        } else {
            return std::pow(ff_dipole(q2, *A), 2.);
        }
    }

    void initCrossSectionCDF(double Trmin, double Trmax, const Target& t, std::size_t n);
    [[deprecated("Use initCrossSectionCDF instead")]]
    void init_integral_of_ff2(double Trmin, double Trmax, const Target& t, std::size_t n) {
        initCrossSectionCDF(Trmin, Trmax, t, n);
    }

    Value ff2Integral(Value q2, const Target& t) const {
        return std::exp(gsl_spline_eval(map_cdf.at(t.name),
                                        std::log(q2),
                                        map_acc.at(t.name)));
    }

    Value inverseCDF(Value xi, const Target& t) const {
        Value Trmax = this->maximumRecoilT(t);
        if (Trmax > Tr_max_init) {
            throw std::runtime_error {
                "The maximum recoil energy exceeds the interval of "
                "interpolation at initialization. Consider increasing the "
                "boundary of interpolation (usually set to the maximal "
                "incident kinetic energy in the simulation)."
            };
        }
        Value cdf = xi * ff2Integral(2*t.m*Trmax, t);
        if (cdf < ff2Integral(2*t.m*Tr_min_init, t)) {
            return Tr_min_init;
        }
        return std::exp(gsl_spline_eval(map_cdf_inv.at(t.name), std::log(cdf),
                                        map_acc_inv.at(t.name))) / (2.0 * t.m);
    }

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override {
        return inverseCDF(rn.uniform_xi(), t);
    }

    Value totalCrossSection(const Target& t) const override {
        Value q2max = 2. * t.m * this->maximumRecoilT(t);
        return this->sigmaEff(t) * ff2Integral(q2max, t) / q2max;
    }
};

template<typename Vector3, typename Value>
void SIHelmDM<Vector3, Value>::initCrossSectionCDF(
    double Trmin, double Trmax, const Target& t, std::size_t n)
{
    Tr_min_init = Trmin;
    Tr_max_init = Trmax;
    gsl_function F;
    F.function = &integrand;
    F.params = (void *) &t.A;
    std::vector<double> q2 = geomspace(2.*t.m*Trmin, 2.*t.m*Trmax*2, n);
    std::vector<double> ff2_integ(n, 0.0);

    double result, error;
    double lower {0.0};

    // cumulatively integrate on grid q2
    for (std::size_t i = 0; i < q2.size(); ++i) {
        gsl_integration_qag(&F, lower, q2[i], 0, 1e-12, 1000, 4, w, &result, &error);
        ff2_integ[i] = i > 0 ? ff2_integ[i-1] + result : result;
        lower = q2[i];
    }

    // logarithmization and make sure ff2_integ strictly increasing
    std::vector<double> logx;
    std::vector<double> logy;
    logx.push_back(std::log(q2[0]));
    logy.push_back(std::log(ff2_integ[0]));
    for (std::size_t i = 1; i < n; ++i) {
        if (std::log(ff2_integ[i]) > std::log(ff2_integ[i-1])) {
            logx.push_back(std::log(q2[i]));
            logy.push_back(std::log(ff2_integ[i]));
        }
    }
    logx.back() = std::log(q2.back());
    logy.back() = std::log(ff2_integ.back());
    std::size_t m = logx.size();

    if (map_cdf.count(t.name) != 0) {
        gsl_spline_free(map_cdf[t.name]);
        gsl_interp_accel_free(map_acc[t.name]);
    }
    map_cdf[t.name] = gsl_spline_alloc(gsl_interp_linear, m);
    gsl_spline_init(map_cdf[t.name], logx.data(), logy.data(), m);
    map_acc[t.name] = gsl_interp_accel_alloc();

    if (map_cdf_inv.count(t.name) != 0) {
        gsl_spline_free(map_cdf_inv[t.name]);
        gsl_interp_accel_free(map_acc_inv[t.name]);
    }
    map_cdf_inv[t.name] = gsl_spline_alloc(gsl_interp_linear, m);
    gsl_spline_init(map_cdf_inv[t.name], logy.data(), logx.data(), m);
    map_acc_inv[t.name] = gsl_interp_accel_alloc();
}

} // namespace darkprop
