/**
 * @file
 * The Sun model and the corresponding utility functions.
 */

#pragma once

#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <spdlog/spdlog.h>
#include "Const.hpp"
#include "MediumBall.hpp"
#include "Vector3d.hpp"

namespace darkprop {

/**
 * A parameterization of solar proton number density. The unit is cm^-3.
 *
 * The argument \p x here is defined as r / RSun, a radius r over solar radius RSun.
 * \p x ranges in [0,1].
 */
template<typename Value=double>
Value _solar_proton_number_density_fitted(Value x)
{
    if (x < 0 || x > 1) { throw std::range_error { "out of range (0, 1)" }; }
    Value x2 = x * x;
    Value result = 2.8026161070632417e23 *
        ((111.41677480208095 / std::exp(22.5 * x2)) +
            ((2836.0753895197354 * x2) / std::exp(100 * x2)) +
            ((113.28444361648015 * x2 * (1 - x2)) / std::exp(11 * x2)) +
            1.3339061580148568 * std::pow(1 - x2, 2) -
            1.103177258895335 * std::pow(1 - x2, 3));

    return result;
}

/**
 * A parameterization of solar helium number density. The unit is cm^-3.
 *
 * The argument \p x here is defined as r / RSun, a radius r over solar radius RSun.
 * \p x ranges in [0,1].
 */
template<typename Value=double>
Value _solar_He4_number_density_fitted(Value x)
{
    if (x < 0 || x > 1) { throw std::range_error { "out of range (0, 1)" }; }
    Value x2 = x * x;
    Value result = 7.054866757596805e22 *
        ((203.17277527058866 / std::exp(161 * x2)) +
            ((8481.39965741552 * x2) / std::exp(81 * x2)) +
            ((20862.42349223336 * x2 * x2) / std::exp(40 * x2)) +
            ((63.04121714425907 * x2 * (1 - x2)) / std::exp(12 * x2)) +
            0.429816558939158 * std::pow(1 - x2, 2) -
            0.3186089761379736 * std::pow(1 - x2, 3));

    return result;
}

/**
 * Solar proton and helium number density.
 */
template<typename Value>
Value solar_number_density(Value r, const Target& target)
{
    Value density = 0;
    if (target.A == 4) {
        density = _solar_He4_number_density_fitted(r / constants::rSun);
    } else if (target.A == 1) {
        density = _solar_proton_number_density_fitted(r / constants::rSun);
    } else {
        throw std::invalid_argument {"Unknown target: " + target.name};
    }
    return density / units::cm3;
}

/**
 * Gives locale temperature of the Sun.
 */
template<typename Value>
Value solar_temperature(Value r)
{
    const Value x = r / constants::rSun;
    if (x < 0 || x > 1) { throw std::range_error { "out of range (0, 1)" }; }
    return (1.3434490188416173
        + x * (0.22887431352379267
        + x * (-42.96961469626326
        + x * (364.56459466435035
        + x * (-3727.355629736401
        + x * (35084.64267333538
        + x * (-225766.42280423897
        + x * (977714.4309763662
        + x * (-2.935978388409898e6
        + x * (6.245094109813822e6
        + x * (-9.481726012898933e6
        + x * (1.0209161475686733e7
        + x * (-7.613368141732147e6
        + x * (3.738448304337283e6
        + x * (-1.087026233486686e6
        + x * 141766.42421330605))))))))))))))) * units::keV;
}


/*
 * Read data from a CSV file and store it in a vector of blocks.
 * Each block is represented as a vector of vectors.
 * In each block, each vector is {y, x, z} with z the normalized number density integral.
 */
template<typename Value=double> std::vector<std::vector<std::vector<Value>>>
_load_density_integral(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file) { throw std::runtime_error("Failed to open the file." + filename); }

    std::vector<std::vector<std::vector<Value>>> blocks;
    std::vector<std::vector<Value>> current_block;
    std::string line;

    while (std::getline(file, line)) {
        if (std::none_of(line.cbegin(), line.cend(),
                         [](unsigned char c) { return std::isdigit(c); })) {
            if (!current_block.empty()) {
                blocks.push_back(current_block);
                current_block.clear();
            }
            continue;
        }

        std::vector<Value> row;
        std::stringstream ss(line);
        std::string cell;

        while (ss >> cell) {
            Value number = std::stod(cell);
            row.push_back(number);
        }

        current_block.push_back(row);
    }

    if (!current_block.empty()) {
        blocks.push_back(current_block);
    }

    return blocks;
}


// Check if a given point is below the line determined by two points.
template<typename Value>
inline bool _below_line(Value x1, Value y1, Value x2, Value y2, Value x, Value y) {
    return y < y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

template<typename Value>
std::vector<Value> _find_triangle_contains_the_point(
    const std::vector<std::vector<Value>>& line1,
    const std::vector<std::vector<Value>>& line2, Value x, Value y)
{
    for (std::size_t i = 0; i < line1.size(); i++) {
        Value x1 = line1[i][0];
        Value y1 = line1[i][1];
        Value z1 = line1[i][2];

        Value x2 = line2[i][0];
        Value y2 = line2[i][1];
        Value z2 = line2[i][2];
        if (_below_line(x1, y1, x2, y2, x, y)) {

            // if i = 0，return the first triangle
            std::size_t tmp = 0;
            if (i == 0) {
                tmp = 2;
            }

            Value x3 = line1[tmp + i - 1][0];
            Value y3 = line1[tmp + i - 1][1];
            Value z3 = line1[tmp + i - 1][2];

            return {x1, y1, z1, x2, y2, z2, x3, y3, z3};
        }

        if (i == line2.size() - 1) {
            // this should be on line1 because y3 on line2 can be all the same
            Value x3 = line1[i - 1][0];
            Value y3 = line1[i - 1][1];
            Value z3 = line1[i - 1][2];
            return {x1, y1, z1, x2, y2, z2, x3, y3, z3};
        }
        x2 = line2[i + 1][0];
        y2 = line2[i + 1][1];
        z2 = line2[i + 1][2];
        if (_below_line(x1, y1, x2, y2, x, y)) {

            Value x3 = line2[i][0];
            Value y3 = line2[i][1];
            Value z3 = line2[i][2];

            return {x1, y1, z1, x2, y2, z2, x3, y3, z3};
        }
    }

    throw std::runtime_error("The triangle is not found.");
}

template<typename Value>
Value triangular_interpolation(const std::vector<Value>& points, Value x, Value y)
{
    Value x1 = points[0];
    Value y1 = points[1];
    Value z1 = points[2];
    Value x2 = points[3];
    Value y2 = points[4];
    Value z2 = points[5];
    Value x3 = points[6];
    Value y3 = points[7];
    Value z3 = points[8];

    Value denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

    Value weight1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator;
    Value weight2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator;
    Value weight3 = 1.0 - weight1 - weight2;

    return weight1 * z1 + weight2 * z2 + weight3 * z3;
}

/**
 * Perform the interpolation. Return x.
 */
template<typename Value>
Value perform_interpolation(Value y, Value xz, const std::vector<Value>& y_list,
        const std::vector<std::vector<std::vector<Value>>>& blocks)
{
    std::size_t i = y == y_list.back() ? y_list.size() - 1 : std::distance(
            y_list.cbegin(), std::upper_bound(y_list.cbegin(), y_list.cend(), y));

    std::vector<Value> neighbor_points = _find_triangle_contains_the_point(
        blocks[i - 1], blocks[i], y, xz);

    return triangular_interpolation(neighbor_points, y, xz);
}


/**
 * A Sun model. Currently, the Sun only works with SolarDM which has constant cross section.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class Sun : public MediumBall<Vector3, Value>
{
private:
    std::vector<std::vector<std::vector<Value>>> blocks;
    std::vector<std::vector<std::vector<Value>>> blocks_inv;  // TODO remove blocks_inv
    std::vector<Value> y_list;
public:
    explicit Sun(Value r=constants::rSun, Value rf=units::cm)
      : MediumBall<Vector3, Value>(r, rf) {
        this->targets = std::vector<Target> {
            Target({"H", 1, 1, constants::mp}),
            Target({"He", 2, 4, 4*constants::mu})
        };
    }
    explicit Sun(const std::string& filename, Value r=constants::rSun, Value rf=units::cm)
      : Sun(r, rf) { init(filename); }
    ~Sun() {}

    void init(const std::string& filename);
    Value densityIntegral(Value y);
    /**
     * The normalized density integral.
     */
    Value densityIntegralNormalized(Value y, Value x) {
        return perform_interpolation(y, x, this->y_list, this->blocks);
    }
    /**
     * The inverse of the normalized density integral.
     */
    Value densityIntegralInverseNormalized(Value y, Value z) {
        return perform_interpolation(y, z, this->y_list, this->blocks_inv);
    }

    Value propagate(Particle<Vector3, Value>& particle, RandomNumber<Value>& rn) override;
    Target sampleTarget(const Particle<Vector3, Value>& particle,
                        RandomNumber<Value>& rn) const override;
};
/**
 * Load the tabular normalized number density integral function.
 */
template<typename Vector3, typename Value>
void Sun<Vector3, Value>::init(const std::string& filename) {
    y_list.clear();
    blocks = _load_density_integral(filename);
    blocks_inv = blocks;
    for (auto& block : blocks_inv) {
        for (auto& row : block) {
            std::swap(row[1], row[2]);
        }
        y_list.push_back(block[0][0]);
    }
}

/**
 * For certain \p y, the maximal value of the integrated number density,
 * i.e. the value of integrated number density at \f$ x=\sqrt{1-y^2} \f$.
 * unit: cm^-3. Same as number density.
 * \p y ranges in [0,1].
 */
template<typename Vector3, typename Value>
Value Sun<Vector3, Value>::densityIntegral(Value y)
{
    if (y < 0 || y > 1) { throw std::runtime_error("y must be in [0, 1]"); }
    Value y2 = y * y;
    Value sqrt_1_y2 = std::sqrt(1.0 - y2);
    Value result = sqrt_1_y2 * (
              1.7350267911907757e24
            + y2 * (-8.292739451071844e24
            + y2 * (1.6914376934952912e25
            + y2 * (-1.7735302150247247e25
            + y2 * (9.223297343969258e24
            + y2 * (-1.8446594687938507e24)))))
        )
        + 1.5995321645518595e25 * std::erf(4.69041575982343 * sqrt_1_y2)
            / std::exp(22.0 * y2)
        + 5.564594520975844e24 * std::erf(8.94427190999916 * sqrt_1_y2)
            / std::exp(80.0 * y2)
        + 1.9558180186014135e24 * std::erf(14.142135623730951 * sqrt_1_y2)
            / std::exp(200.0 * y2);
    return result / units::cm3;
}

template<typename Vector3, typename Value>
Value Sun<Vector3, Value>::propagate(Particle<Vector3, Value>& p, RandomNumber<Value>& rn)
{
    Value x0 = p.r.dot(p.ep) / this->radius;
    Value y = (p.r.cross(p.ep)).norm() / this->radius;
    // The particle can be on the surface only if it points inwards (y != 1)
    if (x0 * x0 + y * y > 1 || y == 1) {
        p.in_medium = false;
        p.r += p.ep * this->rfinal;
        p.updateEr();
        p.t += this->rfinal / p.v.norm();
        spdlog::warn("Particle is out of the Sun with y = {0:.16e}", y);
        return 1.0;
    }
    Value xi = rn.uniform_xi();
    Value density_norm = densityIntegral(y);
    Value z = densityIntegralNormalized(y, x0)
              - log(xi) / (p.sigma0 * this->radius * density_norm);

    Value freep = 0;
    if (z > 1) {
        freep = 2.1 * this->radius;
    } else {
        freep = (densityIntegralInverseNormalized(y, z) - x0) * this->radius;
    }

    if (!std::isfinite(freep)) {
        spdlog::warn("freep = {} is not finite", freep);
        freep = 2.1 * this->radius;
    }

    if (freep < 0) {
        spdlog::warn("freep = {0:.16e} Rsun, where T = {1:.16e} MeV, "
                     "r = ({2:.16e}, {3:.16e}, {4:.16e})",
                     freep / this->radius, p.T / units::MeV, p.r.x / this->radius,
                     p.r.x / this->radius, p.r.x / this->radius);
        freep = std::abs(freep);
    }

    go_straight_and_check_the_boundary(p, freep, this->radius, this->rfinal);
    return 1.0;
}

template<typename Vector3, typename Value>
Target Sun<Vector3, Value>::sampleTarget(
        const Particle<Vector3, Value>& p,
        RandomNumber<Value>& rn) const
{
    Value r = p.r.norm();
    // number density n times cross sections (1 and 4)
    Value ns_p = solar_number_density(r, this->targets[0]);
    Value ns_he4 = solar_number_density(r, this->targets[1]) * 4.0;
    Value ns_sum = ns_p + ns_he4;
    if (rn.uniform_xi() < ns_p / ns_sum) {
        return this->targets[0];
    }
    return this->targets[1];
}

} // namespace darkprop
