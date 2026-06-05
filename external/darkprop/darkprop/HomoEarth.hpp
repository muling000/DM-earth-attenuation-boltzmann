/**
 * @file
 * A homogeneous Earth model.
 */

#pragma once

#include <cmath>
#include <array>
#include <spdlog/spdlog.h>
#include "Const.hpp"
#include "Earth.hpp"
#include "SIDM.hpp"
#include "Vector3d.hpp"

namespace darkprop {

/**
 * A homogeneous Earth model composed of 8 components.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class HomoEarth : public Earth<Vector3, Value>
{
private:
    static constexpr Value density { 2.7 * units::g_cm3 };
    static constexpr std::array<Value, 8> target_density_array {
        46.6e-2 * density / (16*constants::mu),
        27.7e-2 * density / (28*constants::mu),
        8.1e-2  * density / (27*constants::mu),
        5.0e-2  * density / (56*constants::mu),
        3.6e-2  * density / (40*constants::mu),
        2.8e-2  * density / (39*constants::mu),
        2.6e-2  * density / (23*constants::mu),
        2.1e-2  * density / (24*constants::mu)
    };
    bool use_cache=false;
    Value cached_lambda = 0;
    std::array<Value, 8> cached_prob_array;
public:
    HomoEarth() : Earth<Vector3, Value>() {
        this->targets = std::vector<Target> {
            Target({"O", 8, 16, 16*constants::mu}),
            Target({"Si", 14, 28, 28*constants::mu}),
            Target({"Al", 13, 27, 27*constants::mu}),
            Target({"Fe", 26, 56, 56*constants::mu}),
            Target({"Ca", 20, 40, 40*constants::mu}),
            Target({"K", 19, 39, 39*constants::mu}),
            Target({"Na", 11, 23, 23*constants::mu}),
            Target({"Mg", 12, 24, 24*constants::mu})};
    }
    ~HomoEarth() {}
    void setCache(const Particle<Vector3, Value>& p) override;
    inline Value meanFreePath(const Particle<Vector3, Value>& p, std::size_t i) const;
    inline Value meanFreePath(const Particle<Vector3, Value>& p) const;
    Value propagate(Particle<Vector3, Value>& particle,
                    RandomNumber<Value>& rn) override;
    Target sampleTarget(const Particle<Vector3, Value>& particle,
                        RandomNumber<Value>& rn) const override;
};

template<typename Vector3, typename Value>
Value HomoEarth<Vector3, Value>::meanFreePath(
        const Particle<Vector3, Value>& p, size_t i) const
{
    Value lambda = target_density_array[i]*p.totalCrossSection(this->targets[i]);
    return lambda > 0 ? 1.0 / lambda : std::numeric_limits<Value>::infinity();
}

template<typename Vector3, typename Value>
Value HomoEarth<Vector3, Value>::meanFreePath(const Particle<Vector3, Value>& p) const
{
    Value lambda {0.0};
    for (size_t i = 0; i != 8; ++i) {
        lambda += target_density_array[i]*p.totalCrossSection(this->targets[i]);
    }
    return lambda > 0 ? 1.0 / lambda : std::numeric_limits<Value>::infinity();
}

template<typename Vector3, typename Value>
void HomoEarth<Vector3, Value>::setCache(const Particle<Vector3, Value>& p)
{
    cached_lambda = meanFreePath(p);
    for (size_t i = 0; i != 8; ++i) {
        cached_prob_array[i] = cached_lambda / meanFreePath(p, i);
    }
    use_cache = true;
}

template<typename Vector3, typename Value>
Value HomoEarth<Vector3, Value>::propagate(
    Particle<Vector3, Value>& p,
    RandomNumber<Value>& rn)
{
    Value xi = rn.uniform_xi();
    Value freep {0.0};
    Value lambda {use_cache ? cached_lambda : meanFreePath(p)};
    if (std::isfinite(lambda)) {
        freep = -std::log(1 - xi) * lambda;
    } else {
        freep = 2.1 * this->radius;  // leave the Earth
    }
    go_straight_and_check_the_boundary(p, freep, this->radius, this->rfinal);
    return 1.0;
}

template<typename Vector3, typename Value>
Target HomoEarth<Vector3, Value>::sampleTarget(
        const Particle<Vector3, Value>& p,
        RandomNumber<Value>& rn) const
{
    Value xi = rn.uniform_xi();
    Value prob = 0.0;
    if (use_cache) {
        for (size_t i = 0; i != 8; ++i) {
            prob += cached_prob_array[i];
            if (xi < prob) {
                return this->targets[i];
            }
        }
    } else {
        for (size_t i = 0; i != 8; ++i) {
            prob += (meanFreePath(p) / meanFreePath(p, i));
            if (!std::isfinite(prob)) {
                spdlog::warn("prob = {} is not finite", prob);
            }
            if (xi < prob) {
                return this->targets[i];
            }
        }
    }
    // this should not happen
    throw std::runtime_error("No target can be sampled.");
}

} // namespace darkprop
