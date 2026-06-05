/**
 * @file
 * Useful numerical tools.
 */

#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include "Base.hpp"
#include "Const.hpp"

namespace darkprop {

/* Determine the index of a flattened multidimensional array.
 *
 * @param indices The indices of the point in the \p N dimensional array.
 * @param dims Length of each dimension.
 */
template<std::size_t N, typename I, typename D>
constexpr std::size_t flat_index(const I& indices, const D& dims) {
    if constexpr (N == 0) {
        return 0;
    } else {
        return indices[N - 1] + dims[N - 1] * flat_index<N - 1>(indices, dims);
    }
}

/* Generate a Vector of geometric sequence,
 *
 * @param a The starting value of the sequence.
 * @param b The final value of the sequence.
 * @param n Number of samples to generate.
 */
template<typename Value>
constexpr std::vector<Value> geomspace(Value a, Value b, std::size_t n)
{
    if (n <= 1) throw(std::invalid_argument("n <= 1"));
    std::vector<Value> res(n);
    Value q = std::pow(b/a, 1./(n-1));
    res[0] = a;
    for (std::size_t i = 1; i < n; ++i) {
        res[i] = res[i-1] * q;
    }
    return res;
}

/* Generate a Vector of arithmetic sequence.
 *
 * @param a The starting value of the sequence.
 * @param b The final value of the sequence.
 * @param n Number of samples to generate.
 * @param useend If true, \p b is the last sample. Otherwise, it is not included.
 */
template<typename Value>
constexpr std::vector<Value> linspace(Value a, Value b, std::size_t n, bool useend=true)
{
    if (n <= 1) throw(std::invalid_argument("n <= 1"));
    std::vector<Value> res(n);
    Value d = useend ? (b - a) / (n - 1) : (b - a) / n;
    res[0] = a;
    for (std::size_t i = 1; i < n; ++i) {
        res[i] = res[i - 1] + d;
    }
    return res;
}

/**
 * Rotate a vector \p v around another vector \p n by angle \p theta.
 */
template<typename Vector3, typename Value>
void rotation_axis_angle(Vector3& v, const Vector3 n, Value theta)
{
  v = v.dot(n) * n + std::cos(theta) * n.cross(v).cross(n) + std::sin(theta) * n.cross(v);
}


template<typename Vector3, typename Value>
void scatter_fixed_target(Particle<Vector3, Value>& p,
                          Value T,
                          Value costheta,
                          RandomNumber<Value>& rn)
{
    if (costheta * costheta > 1 || std::isnan(costheta)) {
        spdlog::warn("costheta = {0:.18e}, where p.T = {1:.18e}, T' = {2:.18e} ",
                     costheta, p.T, T);
        throw std::runtime_error("costheta > 1 or is nan");
    }
    Vector3 en = {p.ep[2], p.ep[2], -(p.ep[0] + p.ep[1])};
    en.normalize();

    Vector3 n = costheta * p.ep + std::sqrt(1 - costheta*costheta) * en;

    // update T, p3
    p.T = T;
    p.p3 = n * p.pFromT();

    // give a random azimuthal angle 0 ~ 2pi
    Value phi = rn.uniform_phi();
    rotation_axis_angle(p.p3, p.ep, phi);

    p.updateEp();
    p.updateVwithP3T();
}


template<typename Vector3, typename Value>
void go_straight_and_check_the_boundary(Particle<Vector3, Value>& p,
                                        Value freep,
                                        Value rbound,
                                        Value rfinal)
{
    Vector3 rnow = p.r;
    p.r += p.ep * freep;
    if (p.r.norm() > rbound) {
        p.in_medium = false;
        Value rnorm = rnow.norm();
        Value costheta = p.ep.dot(-p.er);  // here p.er has not been updated yet
        freep = rfinal + rnorm * costheta
            + std::sqrt(rnorm*rnorm * (costheta*costheta - 1) + rbound*rbound);
        p.r = rnow + p.ep * freep;
    }
    p.t += freep / p.v.norm();
    p.updateEr();
}

template<typename Value> inline Value reduce_m(Value m1, Value m2)
{
    return m1 * m2 / (m1 + m2);
}

} // namespace darkprop
