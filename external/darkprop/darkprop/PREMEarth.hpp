/**
 * @file
 * A preliminary reference Earth model (PREM). Refer to
 * https://doi.org/10.1016/0031-9201(81)90046-7, arXiv:1312.1202, and arXiv:1706.02249.
 */

#pragma once

#include <stdexcept>
#include <spdlog/spdlog.h>
#include "Earth.hpp"
#include "Const.hpp"
#include "Vector3d.hpp"
#include "numerical/RootFinding.hpp"

namespace darkprop {

/**
 * Preliminary reference Earth model (PREM).
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class PREMEarth : public Earth<Vector3, Value>
{
private:
    const std::array<Value, 10> layer_radii {
        1221.5 * units::km, 3480 * units::km, 5701   * units::km, 5771 * units::km,
        5971   * units::km, 6151 * units::km, 6346.6 * units::km, 6356 * units::km,
        6368   * units::km, constants::rEarth};

    const std::array<Value, 10> a0 {13.0885, 12.5815, 7.9565, 5.3197, 11.2494,
                                    7.1089, 2.6910, 2.9, 2.6, 2.6};
    const std::array<Value, 10> a1 {0, -1.2638, -6.4761, -1.4836, -8.0298,
                                    -3.8045, 0.6924, 0, 0, 0};
    const std::array<Value, 10> a2 {-8.8381, -3.6426, 5.5283, 0, 0,
                                    0, 0, 0, 0, 0};
    const std::array<Value, 10> a3 {0, -5.5281, -3.0807, 0, 0,
                                    0, 0, 0, 0, 0};

    const std::vector<Target> elements_core {
        Target({"Fe", 26, 56, 56*constants::mu}),  // Iron
        Target({"Si", 14, 28, 28*constants::mu}),  // Silicon
        Target({"Ni", 28, 58, 58*constants::mu}),  // Nickel
        Target({"S",  16, 32, 32*constants::mu}),  // Sulfur
        Target({"Cr", 24, 52, 52*constants::mu}),  // Chromium
        Target({"Mn", 25, 55, 55*constants::mu}),  // Manganese
        Target({"P",  15, 31, 31*constants::mu}),  // Phosphorus
        Target({"C",   6, 12, 12*constants::mu}),  // Carbon
        Target({"H",   1,  1,    constants::mp})   // Hydrogen
    };
    static constexpr std::array<Value, 9> elements_core_fraction {
        0.855, 0.06, 0.052, 0.019, 0.009, 0.0003, 0.002, 0.002, 0.0006
    };
    const std::vector<Target> elements_mantle {
        Target({"O",   8, 16, 16*constants::mu}),  // Oxygen
        Target({"Mg", 12, 24, 24*constants::mu}),  // Magnesium
        Target({"Si", 14, 28, 28*constants::mu}),  // Silicon
        Target({"Fe", 26, 56, 56*constants::mu}),  // Iron
        Target({"Ca", 20, 40, 40*constants::mu}),  // Calcium
        Target({"Al", 13, 27, 27*constants::mu}),  // Aluminium
        Target({"Na", 11, 23, 23*constants::mu}),  // Sodium
        Target({"Cr", 24, 52, 52*constants::mu}),  // Chromium
        Target({"Ni", 28, 58, 58*constants::mu}),  // Nickel
        Target({"Mn", 25, 55, 55*constants::mu}),  // Manganese
        Target({"S",  16, 32, 32*constants::mu}),  // Sulfur
        Target({"C",   6, 12, 12*constants::mu}),  // Carbon
        Target({"H",   1,  1,    constants::mp}),  // Hydrogen
        Target({"P",  15, 31, 31*constants::mu})   // Phosphorus
    };
    static constexpr std::array<Value, 14> elements_mantle_fraction {
        0.440, 0.228, 0.21, 0.0626, 0.0253, 0.0235, 0.0027, 0.0026, 0.002,
        0.001, 0.0003, 0.0001, 0.0001, 0.00009
    };
    std::array<Value, 9> scatter_probability_core;
    std::array<Value, 14> scatter_probability_mantle;
    std::array<Value, 2> g_PREM {0, 0};

    Value gFactor(int layer) {
        if (layer <= 1) { return g_PREM[0]; }
        if (layer <= 9) { return g_PREM[1]; }
        return 0.0;
    }

public:
    bool cross_section_energy_dependent;
    Value accuracy;
    std::size_t maxit;
    explicit PREMEarth(bool t_energy_dependent=false,
                       Value t_accuracy=1e-12,
                       std::size_t t_maxit=100)
        : Earth<Vector3, Value>(), cross_section_energy_dependent {t_energy_dependent},
          accuracy {t_accuracy}, maxit {t_maxit} {
        this->targets = elements_mantle;
    }
    ~PREMEarth() {}
    int getLayer(Value r) const {
        return std::lower_bound(layer_radii.begin(), layer_radii.end(), r)
               - layer_radii.begin();
    }
    Value getLayerRadius(int layer) const { return layer_radii[layer]; }
    void setPREM(const Particle<Vector3, Value>& p);
    void setCache(const Particle<Vector3, Value>& p) override { setPREM(p); };
    Value propagate(Particle<Vector3, Value>& particle, RandomNumber<Value>& rn) override;
    Target sampleTarget(const Particle<Vector3, Value>& particle,
                        RandomNumber<Value>& rn) const override;
    Value toNextBoundary(const Vector3& r0, const Vector3& ep,
                         int& on_bound, int& in_layer);
    Value densityIntegral(Value r0, Value cosa, Value L, std::size_t layer);
    Value density(Value r, int layer) {
        Value x = r / this->radius;
        return (a0[layer]+x*(a1[layer]+x*(a2[layer]+x*a3[layer])))*units::g_cm3;
    }
};


template<typename Vector3, typename Value>
void PREMEarth<Vector3, Value>::setPREM(const Particle<Vector3, Value> & p)
{
    g_PREM[0] = 0.0;
    g_PREM[1] = 0.0;

    // core
    for (std::size_t i = 0; i < elements_core.size(); ++i) {
        scatter_probability_core[i] = elements_core_fraction[i]
            / elements_core[i].m * p.totalCrossSection(elements_core[i]);
        g_PREM[0] += scatter_probability_core[i];
    }
    for (Value& prob : scatter_probability_core) { prob /= g_PREM[0]; }

    // mantle
    for (std::size_t i = 0; i < elements_mantle.size(); ++i) {
        scatter_probability_mantle[i] = elements_mantle_fraction[i]
            / elements_mantle[i].m * p.totalCrossSection(elements_mantle[i]);
        g_PREM[1] += scatter_probability_mantle[i];
    }
    for (Value& prob : scatter_probability_mantle) { prob /= g_PREM[1]; }
}

template<typename Vector3, typename Value>
Target PREMEarth<Vector3, Value>::sampleTarget(
        const Particle<Vector3, Value>& p,
        RandomNumber<Value>& rn) const
{
    // Note that setPREM is not called here
    int layer = getLayer(p.r.norm());
    Value xi = rn.uniform_xi();
    Value sum = 0.0;
    // Mantle
    if (layer > 1 && layer < 10) {
        for (std::size_t i = 0; i < elements_mantle.size(); ++i) {
            sum += scatter_probability_mantle[i];
            if (sum > xi) {
                return elements_mantle[i];
            }
        }
    }
    // Core
    else if (layer < 2) {
        for (std::size_t i = 0; i < elements_core.size(); ++i) {
            sum += scatter_probability_core[i];
            if (sum > xi) {
                return elements_core[i];
            }
        }
    }
    throw std::runtime_error("no target");
}

/**
 * Analytic integration of the density.
 */
template<typename Vector3, typename Value>
Value PREMEarth<Vector3, Value>::densityIntegral(Value r0, Value cosa, Value L,
                                                 std::size_t layer)
{
    // TODO make sure f(0) = 0
    if (L == 0) return 0;
    Value r = std::sqrt(L * (L + 2 * r0 * cosa) + r0 * r0);

    Value radius = this->radius;
    Value C0, C1, C2, C3;
    C0 = a0[layer] == 0 ? 0 : L * radius * radius;
    C2 = a2[layer] == 0 ? 0 : L * (r0 * (r0 + L * cosa) + L * L / 3);
    // TODO 1e-14
    if (std::abs(cosa + 1) > 1e-14) {
        Value logfactor = std::log((L + r + r0 * cosa)
                                   / (r0 * (1 + cosa)));
        C1 = a1[layer] == 0 ? 0
            : radius / 2 * (r * (L + r0 * cosa)
                + r0 * r0 * (-cosa + (1 - cosa * cosa) * logfactor));
        C3 = a3[layer] == 0 ? 0
            : 1.0/8.0/radius * ((5 - 3*cosa*cosa) * (r - r0) * r0*r0*r0 * cosa
                + 2 * L * L * r * (L + 3 * r0 * cosa)
                + L * r * r0 * r0 * (5 + cosa * cosa)
                + 3*r0*r0*r0*r0 * std::pow(1 - cosa * cosa, 2) * logfactor);
    } else {
        C1 = radius * (L < r0 ? L * (r0 - L/2.0) : L * (L/2.0 - r0) + r0*r0);
        C3 = 0.25/radius * (r0*r0*r0*r0 - std::abs(r0-L) * std::pow(r0-L, 3));
    }
    return (a0[layer] * C0 + a1[layer] * C1 + a2[layer] * C2 + a3[layer] * C3)
           / radius / radius * units::g_cm3;
}

template<typename Vector3, typename Value>
Value PREMEarth<Vector3, Value>::toNextBoundary(const Vector3& r0,
                                                const Vector3& ep,
                                                int& on_bound,
                                                int& in_layer)
{
    Value r = r0.norm();
    Value cosa = r0.dot(ep) / r;
    Value factor = r * r * (cosa * cosa - 1);
    int layer = getLayer(r);
    Value R;
    // outward
    if (cosa >= 0) {
        on_bound = on_bound < 0 ? layer : on_bound + 1;
    } else if (on_bound != 0 && layer != 0) { // inward
        on_bound = on_bound < 0 ? layer : on_bound;
        R = layer_radii[on_bound - 1];
        Value delta = factor + R * R;
        // inner layer
        if (delta > 0) {
            in_layer = on_bound;
            --on_bound;
            return -r * cosa - std::sqrt(delta);
        }
    } // on_bound == 0 is the same
    in_layer = on_bound;
    R = layer_radii[on_bound];
    return -r * cosa + std::sqrt(factor + R * R);
}

template<typename Vector3, typename Value>
Value PREMEarth<Vector3, Value>::propagate(Particle<Vector3, Value>& p,
                                           RandomNumber<Value>& rn)
{
    if (cross_section_energy_dependent) setPREM(p);
    Value xi = rn.uniform_xi();
    Value nlogxi = -log(1 - xi);
    Vector3 r = p.r;
    int layer = getLayer(r.norm());
    Value r0, cosa, L, intg_new;
    Value intg = 0;
    int on_bound = -1;
    while (on_bound < 9) {
        L = toNextBoundary(r, p.ep, on_bound, layer);
        r0 = r.norm();
        cosa = r.dot(p.ep) / r0;
        intg_new = intg + gFactor(layer) * densityIntegral(r0, cosa, L, layer);
        if (intg_new < nlogxi) {
            // No scattering in this layer
            intg = intg_new;
            r += L * p.ep;
            if (on_bound == 9) {
                r += this->rfinal * p.ep;
                p.in_medium = false;
            }
        } else {
            auto f = [=](Value rts) {
                return intg + gFactor(layer) * densityIntegral(r0, cosa, rts, layer)
                       + std::log(1.0 - xi);
            };
            auto df = [=](Value rts) {
                return gFactor(layer)*density(std::sqrt(rts*(rts+2*r0*cosa)+r0*r0), layer);
            };
            r += numerical::rtsafe<Value>(f, df, 0.0, L, 0.0, accuracy, maxit) * p.ep;
            if (r.norm() > this->radius) {
                // This may happen due to numerical inaccuracy
                r += this->rfinal * p.ep;
                p.in_medium = false;
            }
            break;
        }
    }
    p.t += (r - p.r).norm() / p.v.norm();
    p.r = r;
    p.updateEr();
    return 1.0;
}

} // namespace darkprop
