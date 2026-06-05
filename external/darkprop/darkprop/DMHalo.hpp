/**
 * @file
 * Dark matter halo properties.
 */
#pragma once
#include <stdexcept>
#include <cmath>
#include "Const.hpp"

namespace darkprop {

/**
 * The normalization of Maxwell velocity distribution (vesc = +inf).
 */
template<typename Value=decltype(constants::v0)>
inline Value norm_maxwell(Value v0=constants::v0) {
    return M_PI * std::sqrt(M_PI) * v0 * v0 * v0;
}

/**
 * The normalization of halo DM velocity distribution with a finite escape velocity.
 */
template<typename Value=decltype(constants::vesc)>
inline Value norm_esc(Value vesc=constants::vesc, Value v0=constants::v0) {
    return norm_maxwell(v0) * (std::erf(vesc / v0)
                -M_2_SQRTPI * vesc / v0 * std::exp(-vesc * vesc / v0 / v0));
}

/**
 * Halo DM 1D speed distribution.
 */
template<typename Value>
Value fv_halo_1d(Value v, Value vesc=constants::vesc, Value v0=constants::v0) {
    return 4.0 * M_PI / norm_esc(vesc, v0) * v*v * std::exp(-v * v / v0 / v0);
}

/**
 * Halo DM 1D speed distribution in the earth frame.
 */
template<typename Value>
Value fv_halo_1d_earth(Value v, Value vearth=constants::vearth,
                       Value vesc=constants::vesc, Value v0=constants::v0) {
    if (vearth <= 0) { throw std::invalid_argument("vearth must > 0"); }
    if (v > vesc + vearth) { return 0; }
    Value factor = vesc - vearth > v ? std::exp(-std::pow((v + vearth)/v0, 2))
                                     : std::exp(-std::pow(vesc / v0, 2));
    return M_PI * v0 * v0 * v / (norm_esc(vesc, v0) * vearth)
           * (std::exp(-std::pow((v - vearth)/v0, 2)) - factor);
}

} // namespace darkprop
