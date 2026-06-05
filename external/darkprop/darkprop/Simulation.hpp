/**
 * @file
 * Some useful functions to perform simulations.
 */

#pragma once

#include <spdlog/spdlog.h>
#include "MediumBall.hpp"
#include "Const.hpp"
#include "Analysis.hpp"

namespace darkprop {

/**
 * Run a simulation and collect events that cross a sphere.
 */
template<typename Vector3, typename Value>
std::vector<Event<Vector3, Value>> simulate_cross_sphere(
        Particle<Vector3, Value>& p,
        MediumBall<Vector3, Value>& ball,
        Value Tcut,         ///< cutoff energy
        const Vector3& r0,  ///< center of the sphere
        Value R,            ///< radius
        RandomNumber<Value>& rn,
        Value weight0=1.0,
        unsigned long t_long_track=1000000)
{
    std::vector<Event<Vector3, Value>> events;
    unsigned long scatter_num = 0;
    unsigned long long_track = t_long_track;
    Value weight = weight0;
    auto prev_p = p.toEvent(weight, false);
    do {
        prev_p = p.toEvent(weight, false);
        weight *= propagate(p, ball, rn);
        sphere_cross(r0, R, prev_p, p.toEvent(weight, false), events);
        if (p.in_medium) {
            weight *= scatter(p, ball, rn);
            scatter_num++;
            if (scatter_num > long_track) {
                spdlog::warn("scatter more than {0:d} times; depth = {1:.3f} km",
                             long_track, (ball.getRadius() - p.r.norm())/units::km);
                long_track += t_long_track;
            }
            if (std::isnan(p.r.norm())) {
                spdlog::error("p.r.norm = {0:.16e} km, "
                              "p.v.norm = {1:.16e} km/s, "
                              "p.p3.norm = {2:.16e} GeV",
                              p.r.norm() / units::km,
                              p.v.norm() * units::sec / units::km,
                              p.p3.norm());
                throw std::runtime_error("p.r.norm() is nan");
            }
            if (p.T < Tcut) {
                p.in_medium = false;
            }
        }
    } while (p.in_medium);
    return events;
}


/**
 * Run a simulation and collect crossing events on several depths.
 */
template<typename Vector3, typename Value>
std::vector<std::vector<Event<Vector3, Value>>> simulate_cross_depth(
        Particle<Vector3, Value>& p,
        MediumBall<Vector3, Value>& ball,
        Value Tcut,                ///< cutoff energy
        std::vector<Value> depth,  ///< collect events on these depth
        RandomNumber<Value>& rn,
        Value weight0=1.0,
        unsigned long t_long_track=1000000)
{
    std::vector<std::vector<Event<Vector3, Value>>> events(depth.size());
    unsigned long scatter_num = 0;
    unsigned long long_track = t_long_track;
    Value weight = weight0;
    Vector3 r0(0., 0., 0.);
    auto prev_p = p.toEvent(weight, false);
    do {
        prev_p = p.toEvent(weight, false);
        weight *= propagate(p, ball, rn);
        for (size_t i = 0; i != depth.size(); ++i) {
            sphere_cross(r0, ball.getRadius()-depth[i],
                         prev_p, p.toEvent(weight, false), events[i]);
        }
        if (p.in_medium) {
            weight *= scatter(p, ball, rn);
            scatter_num++;
            if (scatter_num > long_track) {
                spdlog::warn("scatter more than {0:d} times; depth = {1:.3f} km",
                             long_track, (ball.getRadius() - p.r.norm())/units::km);
                long_track += t_long_track;
            }
            if (p.T < Tcut) {
                p.in_medium = false;
            }
        }
    } while (p.in_medium);
    return events;
}


/**
 * Run a simulation and return the complete trajectory.
 */
template<typename Vector3, typename Value>
std::vector<Event<Vector3, Value>> simulate_track(
        Particle<Vector3, Value>& p,
        Medium<Vector3, Value>& m,
        Value Tcut,  ///< cutoff energy
        RandomNumber<Value>& rn,
        Value weight0=1.0,
        unsigned long t_long_track=1000000)
{
    std::vector<Event<Vector3, Value>> events;
    unsigned long scatter_num = 0;
    unsigned long long_track = t_long_track;
    Value weight = weight0;
    do {
        events.push_back(p.toEvent(weight, true));
        weight *= propagate(p, m, rn);
        if (p.in_medium) {
            weight *= scatter(p, m, rn);
            scatter_num++;
            if (scatter_num > long_track) {
                spdlog::warn("scatter more than {0:d} times; r = {1:.3f} km",
                             long_track, p.r.norm() / units::km);
                long_track += t_long_track;
            }
            if (p.T < Tcut) {
                p.in_medium = false;
            }
        }
    } while (p.in_medium);
    events.push_back(p.toEvent(weight, true));
    return events;
}

template<typename Vector3, typename Value>
std::vector<std::size_t> analyse_track(
        Particle<Vector3, Value>& p,
        MediumBall<Vector3, Value>& ball,
        Value Tcut,
        Value depth,
        RandomNumber<Value>& rn)
{
    std::vector<Event<Vector3, Value>> events;
    std::size_t scatter_num = 0;
    std::size_t free_particles = 0;
    std::size_t below_vcut = 0;
    Value weight = 1.0;
    do {
        auto prev_p = p.toEvent(weight, false);
        weight *= propagate(p, ball, rn);
        sphere_cross({0, 0, 0}, ball.getRadius() - depth,
                     prev_p, p.toEvent(weight, false), events);

        if (p.in_medium) {
            weight *= scatter(p, ball, rn);
            scatter_num++;
            if (p.T < Tcut) {
                p.in_medium = false;
                below_vcut++;
            }
        }
    } while (p.in_medium);
    std::size_t depth_crossing = events.size();
    if (scatter_num == 0) {
        free_particles++;
    }
    return std::vector<std::size_t> {below_vcut, free_particles,
                                     scatter_num, depth_crossing};
}

} // namespace darkprop
