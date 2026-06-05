/**
 * @file
 * Simple elastic scattering dark matter with Helm form factor.
 */

#include <map>
#include <gsl/gsl_integration.h>
#include "SIDM.hpp"
#include "Utils.hpp"
#include "FormFactor.hpp"

namespace darkprop {

/**
 * Spin independent elastic scattering DM model with nucleus form factor.
 * No interpolation version.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class SIHelmDiscreteDM : public SIDM<Vector3, Value>
{
private:
    gsl_integration_workspace *w;
    std::map<std::string, std::vector<Value>> map_integral_of_ff2;
    std::map<std::string, std::vector<Value>> map_array_q2;

public:
    SIHelmDiscreteDM() : SIDM<Vector3, Value>() {
        w = gsl_integration_workspace_alloc(1000);
    }

    ~SIHelmDiscreteDM() {
        gsl_integration_workspace_free(w);
    }

    static double integrand(double q2, void *t_A) {
        int *A = static_cast<int *>(t_A);
        return std::pow(ff_Helm(q2, *A), 2.);
    }

    static double integrand_log(double logq2, void *t_A) {
        double q2 = std::exp(logq2);
        return integrand(q2, t_A) * q2;
    }

    void init_integral_of_ff2(double Trmin, double Trmax,
                              const Target& t, int n) {
        gsl_function F;
        // F.function = &integrand;
        F.function = &integrand_log;
        F.params = (void *) &t.A;
        std::vector<double> x = geomspace(2.*t.m*Trmin, 2.*t.m*Trmax*2, n);
        std::vector<double> y(n);
        double error;

        // lower than q^2min to avoid zero integral region
        double lower = std::log(x[0] / 10.);

        // log integrate
        for (size_t i = 0; i < x.size(); ++i) {
            gsl_integration_qag(&F, lower, std::log(x[i]), 0, 1e-10, 1000,
                                4, w, &y[i], &error);
        }
        // for (size_t i = 0; i < x.size(); ++i) {
        //     gsl_integration_qag(&F, 0.0, x[i], 0, 1e-6, 1000, 1, w,
        //                         &y[i], &error);
        // }

        // TODO the sort problem
        map_array_q2[t.name] = std::move(x);
        map_integral_of_ff2[t.name] = std::move(y);
    }

    std::size_t _find_q2max_idx(const Target& t) const {
        Value q2max = 2.0 * t.m * this->maximumRecoilT(t);
        auto itqmax = std::lower_bound(map_array_q2.at(t.name).begin(),
                                       map_array_q2.at(t.name).end(), q2max);
        return itqmax - map_array_q2.at(t.name).begin();
    }

    Value scatteringSampleTr(const Target& t, RandomNumber<Value>& rn) override
    {
        Value Tr_max = this->maximumRecoilT(t);
        Value xi = rn.uniform_xi();
        if (t.m * Tr_max < 0.004) {
            // neglect form factor
            // return Tr_max * xi;
            return this->T - Tr_max * xi;
        }
        std::size_t idxqmax = _find_q2max_idx(t);
        Value y = map_integral_of_ff2.at(t.name)[idxqmax] * xi;
        auto it1 = std::lower_bound(map_integral_of_ff2.at(t.name).begin(),
                                    map_integral_of_ff2.at(t.name).end(), y);
        std::size_t idx1 = it1 - map_integral_of_ff2.at(t.name).begin();
        std::size_t idx0 = idx1 - 1;
        return map_array_q2.at(t.name)[idx0] / (2.0 * t.m);
    }

    Value totalCrossSection(const Target& t) const override {
        std::size_t idx = _find_q2max_idx(t);
        Value Trmax = this->maximumRecoilT(t);
        return this->sigmaEff(t) * map_integral_of_ff2.at(t.name)[idx] / (2.*t.m*Trmax);
    }
};

} // namespace darkprop
