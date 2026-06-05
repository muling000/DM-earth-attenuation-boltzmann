#pragma once

#include <vector>
#include <array>
#include <bitset>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <type_traits>
#include "../Utils.hpp"

namespace darkprop::numerical {

/**
 * Linear interpolation on N dimensional regular grid.
 */
template<std::size_t N, typename V>
class GridLinearInterpolator
{
    using T = typename V::value_type;
private:
    V m_f;
    T m_fill_value;
    std::array<std::size_t, N> m_nx;
    std::array<V, N> m_x;
    std::array<T, N> m_x_min;
    std::array<T, N> m_x_max;
public:
    explicit GridLinearInterpolator(T t_fill_value=0.0) : m_fill_value {t_fill_value} {}
    GridLinearInterpolator(const std::array<V, N>& t_x, const V& t_f, T t_fill_value=0.0) {
        interpolate(t_x, t_f, t_fill_value);
    }
    ~GridLinearInterpolator() {}

    void interpolate(const std::array<V, N>& t_x, const V& t_f, T t_fill_value=0.0);

    T getXmin(std::size_t i) const { return m_x_min[i]; }
    T getXmax(std::size_t i) const { return m_x_max[i]; }
    const V& getX(std::size_t i) const { return m_x[i]; }
    const V& getF(std::size_t i) const { return m_f[i]; }
    const auto& getX() const { return m_x; }
    const auto& getF() const { return m_f; }

    T interp(const V& x) const;
    T operator()(const V& x) const { return interp(x); };
    T operator()(const T* x) const { return interp(V(x, x + N)); };
    T operator()(T* x) const { return operator()(const_cast<T*>(x)); };
    // for gsl integrator
    T operator()(double *x, std::size_t dim) const {
        assert(dim == N);
        return operator()(x);
    }
    template<typename... X>
    T operator()(X... x) const {
        static_assert(sizeof...(X) == N);
        return interp({x...});
    }
};

template<std::size_t N, typename V>
GridLinearInterpolator(const std::array<V, N>&, const V&) -> GridLinearInterpolator<N, V>;

template<std::size_t N, typename V>
void GridLinearInterpolator<N, V>::interpolate(
    const std::array<V, N>& t_x, const V& t_f, T t_fill_value)
{
    m_x = t_x;
    m_f = t_f;
    m_fill_value = t_fill_value;
    for (std::size_t i = 0; i < N; ++i) {
        m_nx[i] = m_x[i].size();
        m_x_min[i] = *std::min_element(m_x[i].cbegin(), m_x[i].cend());
        m_x_max[i] = *std::max_element(m_x[i].cbegin(), m_x[i].cend());
    }
}

template<std::size_t N, typename V>
typename GridLinearInterpolator<N, V>::T
GridLinearInterpolator<N, V>::interp(const V& x) const
{
    for (std::size_t i = 0; i < N; ++i) {
        if (x[i] < m_x_min[i] || x[i] > m_x_max[i]) {
            return m_fill_value;
        };
    }
    std::array<std::size_t, N> ix;
    std::array<T, N> r;
    for (std::size_t i = 0; i < N; ++i) {
        ix[i] = x[i] == m_x_max[i] ? m_nx[i] - 1 : std::distance(
            m_x[i].cbegin(), std::upper_bound(m_x[i].cbegin(), m_x[i].cend(), x[i]));
        r[i] = (x[i] - m_x[i][ix[i] - 1]) / (m_x[i][ix[i]] - m_x[i][ix[i] - 1]);
    }
    T value = 0;
    T tmp = 0;
    std::array<std::size_t, N> selected_indices;
    // Loop combinations represented as binary number 0 ~ 2^N.
    // For each bit, 0 represents (ix - 1), 1 represents ix.
    for (std::size_t c = 0; c < 1 << N; ++c) {
        std::bitset<N> bits(c);
        for (std::size_t i = 0; i < N; ++i) {
            selected_indices[i] = bits[i] ? ix[i] : ix[i] - 1;
        }
        tmp = m_f[flat_index<N>(selected_indices, m_nx)];
        for (std::size_t i = 0; i < N; ++i) {
            tmp *= (bits[i] ? r[i] : 1 - r[i]);
        }
        value += tmp;
    }
    return value;
}

} // namespace darkprop::numerical
