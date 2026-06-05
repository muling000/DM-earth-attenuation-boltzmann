#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>

namespace darkprop::numerical {

template<typename F>
class Integrator
{
private:
    gsl_integration_workspace *w;
    gsl_function gF;
public:
    bool use_log;
    int limit, key;
    double epsabs, epsrel, result, error;
    void setLog(bool t_use_log=true) {
        use_log = t_use_log;
        if (use_log) {
            gF.function = &functor_gsl_wrapper_log;
        } else {
            gF.function = &functor_gsl_wrapper;
        }
    }
    static double functor_gsl_wrapper_log(double logx, void* functor)
    {
        F* f = reinterpret_cast<F *>(functor);
        double x = std::exp(logx);
        return (*f)(x) * x;
    }
    static double functor_gsl_wrapper(double x, void* functor)
    {
        F* f = reinterpret_cast<F *>(functor);
        return (*f)(x);
    }
    explicit Integrator(bool t_log=false, double t_epsabs=0.0, double t_epsrel=1e-4,
                        int t_limit=100000, int t_key=4)
        : use_log {t_log}, limit {t_limit}, key {t_key},
          epsabs {t_epsabs}, epsrel {t_epsrel}, result {0}, error {0} {
        w = gsl_integration_workspace_alloc(t_limit);
        setLog(use_log);
    }
    double operator()(const F& f, double t_lo, double t_up) {
        gF.params = (void *) &f;

        double lo = use_log ? std::log(t_lo) : t_lo;
        double up = use_log ? std::log(t_up) : t_up;
        gsl_integration_qag(&gF, lo, up, epsabs, epsrel, limit, key, w, &result, &error);
        return result;
    }
    ~Integrator() {
        gsl_integration_workspace_free(w);
    }
};

/**
 * Wrapper of the GSL multidimensional Monte Carlo integration.
 * The signature of the integrand has to be (double *, std::size_t) -> Value.
 */
template<typename F>
class MCIntegrator
{
private:
    const gsl_rng_type *T;
    gsl_rng *rng;
    gsl_monte_function gF;
public:
    std::size_t calls;
    double result, error;
    explicit MCIntegrator(std::size_t t_calls=100000)
        : calls {t_calls}, result {0}, error {0} {
        gF.f = &functor_wrapper;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        rng = gsl_rng_alloc(T);
    }
    ~MCIntegrator() { gsl_rng_free(rng); }
    static double functor_wrapper(double *x, std::size_t dim, void* functor) {
        F* f = reinterpret_cast<F *>(functor);
        return (*f)(x, dim);
    }
    double operator()(const F& f, std::vector<double> xl, std::vector<double> xu) {
        auto dim = xl.size();
        if (dim != xu.size()) {
            throw std::invalid_argument("xl and xl dimension don't match");
        }
        gF.params = (void *) &f;
        gF.dim = dim;

        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
        do {
            gsl_monte_vegas_integrate(&gF, xl.data(), xu.data(), dim,
                                      calls, rng, s, &result, &error);
        } while (std::abs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);
        gsl_monte_vegas_free(s);
        return result;
    }
};

} // namespace darkprop::numerical
