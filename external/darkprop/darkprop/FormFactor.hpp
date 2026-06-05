/**
 * @file
 * The nucleus form factor.
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include "Const.hpp"

namespace darkprop {

/**
 * This function is used to handle x -> 0.
 */
template<typename Value>
Value _ff_Helm_core(Value x)
{
    if (x >= 1e-1) {
        return 3.*(std::sin(x) - x * std::cos(x)) / (x*x*x);
    }
    if (x > 1e-2) {
        Value x2 = x * x;
        return 1.0 + x2 * (
                -0.1 + x2 * (
                 3.5714285714285714286e-3 + x2 * (
                  -6.6137566137566137566e-5 + x2 * 7.5156325156325156325e-7)));
    }
    if (x > 1e-4) {
        Value x2 = x * x;
        return 1.0 + x2 * (-0.1 + x2 * 3.5714285714285714e-3);
    }
    if (x > 1e-8) {
        return 1.0 - 0.1 * x * x;
    }
    return 1.0;
}


/**
 * Helm form factor.
 *
 * @param qsquare square of momentum transfer
 * @param A mass number
 */
template<typename Value>
Value ff_Helm(Value qsquare, int A)
{
    if (A < 7) {
        throw std::runtime_error("Helm form factor A must be greater than 6");
    }
    Value q = std::sqrt(qsquare);
    Value qr1 = q * std::sqrt(std::pow(1.2 * std::pow(A, 1. / 3.), 2.) - 5.)
        * units::fm;
    return _ff_Helm_core(qr1) * std::exp(-0.5*qsquare * units::fm * units::fm);
    // return 3. * std::sph_bessel(1, qr1) / qr1
    //     * std::exp(-0.5*qsquare * units::fm * units::fm);
}


// template<typename Value>
// Value ff_Helm(Value qsquare, int A) // JCAP04(2007)012
// {
//     Value q = std::sqrt(qsquare);
//     Value qr1 = q * std::sqrt(std::pow(0.91*std::pow(A, 1./3.)+0.3, 2.)-5.)
//         * units::fm;
//     return 3. * std::sph_bessel(1, qr1) / qr1
//         * std::exp(-0.5*qsquare * units::fm * units::fm);
// }

/**
 * Dipole form factor.
 *
 * @param qsquare square of momentum transfer
 * @param A mass number
 */
template<typename Value>
Value ff_dipole(Value qsquare, int A)
{
    Value lambda;
    if (A == 1) {
        lambda = 0.77;
    } else if (A == 4) {
        lambda = 0.41;
    } else {
        throw std::runtime_error("dipole form factor A only supports 1 and 4");
    }
    return 1. / std::pow((1. + qsquare / (lambda * lambda)), 2.0);
}

} // namespace darkprop
