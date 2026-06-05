#pragma once

#include <cmath>
#include <stdexcept>
#include <spdlog/spdlog.h>

namespace darkprop::numerical {

/**
 * The rtsafe algorithm adopted from the book Numerical Recipes.
 *
 * Press, William H., and Saul A. Teukolsky. Numerical recipes 3rd edition:
 * The art of scientific computing. Cambridge university press, 2007.
 */
template<typename Value=double, typename F, typename DF>
Value rtsafe(const F& t_f, const DF& t_df, Value x1, Value x2,
             Value abstol=0.0, Value reltol=1e-12, std::size_t MAXIT=100)
{
    Value fl = t_f(x1);
    Value fh = t_f(x2);
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        spdlog::error("fl = {0}, fh = {1}", fl, fh);
        throw std::runtime_error("Root must be bracketed in rtsafe.");
    }
    if (fl == 0.0) { return x1; }
    if (fh == 0.0) { return x2; }

    Value xl, xh;
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    Value rts = 0.5 * (x1 + x2);
    Value dxold = std::abs(x2 - x1);
    Value dx = dxold;
    Value f = t_f(rts);
    Value df = t_df(rts);
    for (std::size_t i = 0; i < MAXIT; ++i) {
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
            || (std::abs(2.0 * f) > std::abs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) { return rts; } // Change in root is negligible
        } else {
            dxold = dx;
            dx = f / df;
            Value tmp = rts;
            rts -= dx;
            if (tmp == rts) { return rts; }
        }
        if (std::abs(dx) < abstol + reltol * std::abs(rts)) { return rts; }
        f = t_f(rts);
        df = t_df(rts);
        if (f < 0.0) {
            xl = rts;
        } else {
            xh = rts;
        }
    }
    throw std::runtime_error("Number of iterations exceeded in rtsafe");
}

} // namespace darkprop::numerical
