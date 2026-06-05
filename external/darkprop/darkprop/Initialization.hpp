/**
 * @file
 * Initial condition of dark matter.
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <spdlog/spdlog.h>
#include "Base.hpp"
#include "DMHalo.hpp"

namespace darkprop {

/**
 * Sample 1D halo speed from the isotropic SHM.
 */
template<typename Value>
Value sample_halo_speed(RandomNumber<Value>& rn, /**< random number gen */
                        Value vmin=0.0, /**< lower bound */
                        Value vmax=constants::vesc, /**< upper bound */
                        Value vesc=constants::vesc, /**< escape speed */
                        Value v0=constants::v0 /**< most probable speed */)
{
    // TODO use inverse CDF sampling
	Value v, fv;
    do {
        v = rn.uniform_ab(vmin, vmax);
        fv = rn.uniform_xi() * fv_halo_1d(v0, vesc, v0);
    } while (fv > fv_halo_1d(v, vesc, v0));
	return v;
}

/**
 * Sample a random vector isotropically with length \p len (default 1).
 */
template<typename Vector3, typename Value>
Vector3 random_isotropic_vector(RandomNumber<Value>& rn, Value len=1)
{
    Value theta = std::acos(rn.uniform_costheta());
    Value phi = rn.uniform_phi();
    return Vector3(len * std::sin(theta) * std::cos(phi),
                   len * std::sin(theta) * std::sin(phi),
                   len * std::cos(theta));
}

/**
 * Sample 3D halo velocity in the Earth frame from the SHM.
 *
 * The sampling range can be set to (\p vmin, \p vmax). If \p vmax <= 0, it
 * will be set to \p v_earth.norm() + \p vesc.
 */
template<typename Vector3, typename Value>
Vector3 sample_halo_velocity_earth_frame(
        Vector3 v_earth,            /**< earth velocity */
        RandomNumber<Value>& rn,    /**< random number generator */
        Value vmin=0.0,             /**< speed lower bound */
        Value vmax=-1.0,            /**< speed upper bound */
        Value vesc=constants::vesc, /**< halo escape speed */
        Value v0=constants::v0      /**< halo most probable speed */)
{
    Vector3 velocity;
    if (vmax <= 0) {
        vmax = v_earth.norm() + vesc;
    }
    Value halo_vmin = std::max(vmin - v_earth.norm(), static_cast<Value>(0.0));
    Value halo_vmax = std::min(vmax + v_earth.norm(), vesc);
    do {
        Value vi = sample_halo_speed(rn, halo_vmin, halo_vmax, vesc, v0);
        velocity = random_isotropic_vector<Vector3>(rn, vi) - v_earth;
    } while (velocity.norm() < vmin || velocity.norm() > vmax);
    return velocity;
}

/**
 * Sample an incident position according to the given direction.
 */
template<typename Vector3, typename Value>
Vector3 random_incident_position(RandomNumber<Value>& rn,
                                 const Vector3& vi,
                                 Value r_earth=constants::rEarth,
                                 Value R=1.1*constants::rEarth)
{
    Vector3 ez = -vi.normalized();
    Value n = std::sqrt(ez[1] * ez[1] + ez[2] * ez[2]);
    Vector3 ex(0, ez[2] / n, -ez[1] / n);
    Vector3 ey = ez.cross(ex);
    Value phi = rn.uniform_phi();
    Value xi = rn.uniform_xi();
    return R * ez + std::sqrt(xi) * r_earth * (std::cos(phi) * ex + std::sin(phi) * ey);
}


/**
 * Set halo DM random initial condition.
 *
 * \p vmin and \p vmax are in the Earth frame.
 */
template<typename Vector3, typename Value>
void init_halo(Particle<Vector3, Value>& p,
               Value t,                        /**< initial time */
               Value vmin,                     /**< speed lower bound */
               const Vector3& v_earth,         /**< velocity of the Earth */
               RandomNumber<Value>& rn,        /**< random number generator */
               Value R,                        /**< initial disk distance */
               Value vmax=-1.0,                /**< speed upper bound */
               Value vesc=constants::vesc,     /**< halo escape velocity */
               Value v0=constants::v0,         /**< halo most probable speed */
               Value radius=constants::rEarth)
{
    Vector3 vinit = sample_halo_velocity_earth_frame(v_earth, rn, vmin, vmax, vesc, v0);
    Vector3 xinit = random_incident_position(rn, vinit, radius, R);
    p.t = t;
    p.setV(vinit);
    p.setR(xinit);
}


template<typename Value>
inline void check_non_positive(Value T)
{
    if (T <= 0) { throw std::invalid_argument {"parameter <= 0"}; }
}

/**
 * Set the initial condition of a particle with energy T and direction in the galactic
 * coordinate (b, l).
 */
template<typename Vector3, typename Value>
void init_Tbl(Particle<Vector3, Value>& p,
              Value T, Value b, Value l,
              RandomNumber<Value>& rn,
              Value t=0.0,
              Value R=1.1*constants::rEarth,
              Value radius=constants::rEarth)
{
    check_non_positive(T);
    Value pi = std::sqrt(T * (T + 2. * p.m));
    Vector3 ez = Vector3(std::cos(b) * std::cos(l),
                         std::cos(b) * std::sin(l),
                         std::sin(b));
    Vector3 pinit = -pi * ez;
    Vector3 xinit = random_incident_position(rn, pinit, radius, R);
    p.t = t;
    p.setP(pinit);
    p.setR(xinit);
}

template<typename Vector3, typename Value>
void init_Tbl_vertical(Particle<Vector3, Value>& p,
                       Value T, Value b, Value l,
                       Value t=0.0,
                       Value R=constants::rEarth)
{
    check_non_positive(T);
    Value pi = std::sqrt(T * (T + 2. * p.m));
    Vector3 er = Vector3(std::cos(b) * std::cos(l),
                         std::cos(b) * std::sin(l),
                         std::sin(b));
    p.setP(-pi * er);
    p.setR(R * er);
    p.t = t;
}

template<typename Vector3, typename Value>
void init_T_point_source(Particle<Vector3, Value>& p,
                         Value T,
                         const Vector3& ep,
                         RandomNumber<Value>& rn,
                         Value t=0.0,
                         Value R=1.1*constants::rEarth,
                         Value radius=constants::rEarth)
{
    check_non_positive(T);
    Value pi = std::sqrt(T * (T + 2. * p.m));
    Vector3 pinit = pi * ep.normalized();
    Vector3 xinit = random_incident_position(rn, pinit, radius, R);
    p.t = t;
    p.setP(pinit);
    p.setR(xinit);
}

/**
 * Set the initial condition of a particle isotropically with energy T.
 */
template<typename Vector3, typename Value>
void init_T_isotropic(Particle<Vector3, Value>& p,
                      Value T,
                      RandomNumber<Value>& rn,
                      Value t=0.0,
                      Value R=1.1*constants::rEarth,
                      Value radius=constants::rEarth)
{
    check_non_positive(T);
    Value pi = std::sqrt(T * (T + 2. * p.m));
    Vector3 pinit = random_isotropic_vector<Vector3>(rn, pi);
    Vector3 xinit = random_incident_position(rn, pinit, radius, R);
    p.t = t;
    p.setP(pinit);
    p.updateVwithP3T();
    p.setR(xinit);
}

/**
 * Set the initial condition of a particle isotropically on a surface with energy T.
 */
template<typename Vector3, typename Value>
void init_rT_isotropic(Particle<Vector3, Value>& p,
                       Value r,
                       Value T,
                       RandomNumber<Value>& rn,
                       Value t=0.0)
{
    check_non_positive(r);
    check_non_positive(T);

    Value pi = std::sqrt(T * (T + 2. * p.m));
    Vector3 pinit = random_isotropic_vector<Vector3>(rn, pi);
    Vector3 xinit = random_isotropic_vector<Vector3>(rn, r);
    p.setP(pinit);
    p.setR(xinit);
    p.t = t;
    p.in_medium = true;
}

/**
 * Place the particle on a surface.
 */
template<typename Vector3, typename Value>
void inject(Particle<Vector3, Value>& p, Value radius=constants::rEarth)
{
    // already in the medium
    if (p.r.norm() <= radius) {
        p.in_medium = true;
        return;
    }

    Value rep = p.r.dot(p.ep);
    if (rep >= 0) {
        throw std::runtime_error("Particle did not point to the medium.");
    }

    Value delta = rep * rep - p.r.dot(p.r) + radius * radius;

    if (delta <= 0) { throw std::runtime_error("Particle missed the medium."); }

    p.r += (-rep - std::sqrt(delta)) * p.ep;

    if (p.r.norm() > radius) {
        Value x = std::sqrt(delta);
        p.r += (x > 1.0 * units::mm ? 1.0 * units::mm * p.ep : x * p.ep);
    }
    if (p.r.norm() > radius) {
        spdlog::warn("particle injection failed");
        p.updateEr();
        p.in_medium = false;
        return;
    }

    p.updateEr();
    p.in_medium = true;
}

} // namespace darkprop
