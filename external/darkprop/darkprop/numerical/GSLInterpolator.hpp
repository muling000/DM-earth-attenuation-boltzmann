#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <gsl/gsl_spline.h>

namespace darkprop::numerical {

template<typename Value>
class GSLInterpolator
{
private:
    gsl_spline *spline = NULL;
    gsl_interp_accel *acc = NULL;
    std::string method_gsl;
    std::string method;

public:
    explicit GSLInterpolator(const std::string& t_method="linear")
        : method_gsl {t_method} {}
    ~GSLInterpolator();
    void interpolate(const std::vector<Value>& x, const std::vector<Value>& y,
                     const std::string& method="loglog");
    Value operator()(Value);
};

template<typename Value>
GSLInterpolator<Value>::~GSLInterpolator()
{
    if (spline) { gsl_spline_free(spline); }
    if (acc) { gsl_interp_accel_free(acc); }
}

template<typename Value>
void GSLInterpolator<Value>::interpolate(const std::vector<Value>& x,
                                         const std::vector<Value>& y,
                                         const std::string& t_method)
{
    std::size_t N = x.size();
    if (N != y.size()) {
        throw std::runtime_error("x y dimension not match");
    }
    method = t_method;

    // double for gsl
    std::vector<double> interpx(N);
    std::vector<double> interpy(N);

    if (spline) {
        gsl_spline_free(spline);
    }
    spline = gsl_spline_alloc(gsl_interp_linear, N);
    if (acc) {
        gsl_interp_accel_free(acc);
    }
    acc = gsl_interp_accel_alloc();

    // TODO check x == 0 && static_cast
    for (std::size_t i = 0; i != x.size(); ++i) {
        if (method == "loglog") {
            interpx[i] = std::log(x[i]);
            interpy[i] = y[i] > 0 ? std::log(y[i]) : std::log(1e-300);
        } else if (method == "xlog") {
            interpx[i] = std::log(x[i]);
            interpy[i] = y[i];
        } else if (method == "ylog") {
            interpx[i] = x[i];
            interpy[i] = y[i] > 0 ? std::log(y[i]) : std::log(1e-300);
        } else if (method == "linear") {
            interpx[i] = x[i];
            interpy[i] = y[i];
        }
    }
    gsl_spline_init(spline, interpx.data(), interpy.data(), N);
}
template<typename Value>
Value GSLInterpolator<Value>::operator()(Value t_x)
{
    double x = static_cast<double>(t_x);
    double result = 0;
    if (method == "loglog") {
        result = std::exp(gsl_spline_eval(spline, std::log(x), acc));
    } else if (method == "xlog") {
        result = gsl_spline_eval(spline, std::log(x), acc);
    } else if (method == "ylog") {
        result = std::exp(gsl_spline_eval(spline, x, acc));
    } else if (method == "linear") {
        result = gsl_spline_eval(spline, x, acc);
    }
    return static_cast<Value>(result);
}

} // namespace darkprop::numerical
