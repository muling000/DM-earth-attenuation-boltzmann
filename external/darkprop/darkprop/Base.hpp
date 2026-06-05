/**
 * @file
 * Core classes definition.
 */

#pragma once

#include <ios>
#include <cmath>
#include <random>
#include <string>
#include "Const.hpp"

namespace darkprop {

/**
 * Target nucleus or electron.
 */
struct Target
{
    std::string name;  ///< name string
    int Z;             ///< charge number
    int A;             ///< mass number
    double m;          ///< mass
    Target() {}
    Target(const std::string& t_name, int t_Z, int t_A, double t_m)
        : name {t_name}, Z {t_Z}, A {t_A}, m {t_m} {}
    ~Target() {}
};


/**
 * Random number generator.
 */
template<typename Value=double>
class RandomNumber
{
private:
    std::uniform_real_distribution<Value> distrphi;
    std::uniform_real_distribution<Value> distrcos;
    std::uniform_real_distribution<Value> distr01;
    std::mt19937_64 prng;
public:
    explicit RandomNumber(int seed=1) :
        distrphi {std::uniform_real_distribution<Value>(0.0, 2*M_PI)},
        distrcos {std::uniform_real_distribution<Value>(-1., 1.)},
        distr01 {std::uniform_real_distribution<Value>(0., 1.)} {
        if (seed < 0) {
            std::random_device rd;
            prng = std::mt19937_64(rd());
        } else {
            prng = std::mt19937_64(seed);
        }
    }
    /**
     * Sample uniform distribution on (0, pi).
     */
    Value uniform_phi() {
        return distrphi(prng);
    }
    /**
     * Sample uniform distribution on (-1, 1).
     */
    Value uniform_costheta() {
        return distrcos(prng);
    }
    /**
     * Sample uniform distribution on (0, 1).
     */
    Value uniform_xi() {
        return distr01(prng);
    }
    /**
     * Sample uniform distribution on (a, b).
     */
    Value uniform_ab(Value a, Value b) {
        return a + (b - a) * uniform_xi();
    }
    /**
     * Sample normal distribution N(mu, sigma).
     */
    Value normal(Value mu, Value sigma) {
        return std::normal_distribution<Value>(mu, sigma)(prng);
    }
};


/**
 * Store particle information for dump.
 */
template<typename Vector3, typename Value>
struct Event
{
    Value t;       ///< time
    Value T;       ///< kinetic energy
    Vector3 r;     ///< position
    Vector3 p3;    ///< 3 momentum
    Value weight;  ///< statistical weight for importance sampling
    Event() {}
    Event(Value t_t, Value t_T, Vector3 t_r, Vector3 t_p3, Value t_weight)
        : t {t_t}, T{t_T}, r{t_r}, p3{t_p3}, weight{t_weight} {}
    ~Event() {}
};


/**
 * Base abstract Particle class.
 */
template<typename Vector3, typename Value>
class Particle
{
public:
    // TODO make them protected members
    Value m;        ///< mass
    Value sigma0;   ///< cross section
    Value t;        ///< time
    Value T;        ///< kinetic energy
    Vector3 r;      ///< position
    Vector3 v;      ///< velocity
    Vector3 p3;     ///< 3 momentum
    Vector3 ep;     ///< velocity (momentum) direction
    Vector3 er;     ///< position direction
    bool in_medium; ///< switch indicating whether to continue the simulation
public:
    explicit Particle(bool t_in_medium=false) : in_medium {t_in_medium}  {};
    virtual ~Particle() {};

    Value pFromT() const { return std::sqrt(T * (T + 2.0 * m)); }
    Value TFromP3() const {
        Value pnorm = p3.norm();
        if (pnorm / m < 1e-4) {
            return 0.5 * pnorm * pnorm / m;
        }
        return std::sqrt(pnorm * pnorm + m * m) - m;
    }
    void updateEr() {
        er = r.normalized();
    }
    void updateEp() {
        ep = p3.normalized();
    }
    void updateVwithP3T() {
        v = p3 / (T + m);
    }
    void updateP3TwithV() {
        // TODO
        Value vnorm = v.norm();
        if (vnorm < 1e-4) {
            T = 0.5 * m * vnorm * vnorm;
            p3 = m * v;
        } else {
            Value gamma = 1.0 / std::sqrt((1 - vnorm) * (1 + vnorm));
            T = (gamma - 1) * m;
            p3 = gamma * m * v;
        }
        updateEp();
    }
    void setR(Vector3 t_r) {
        this->r = t_r;
        this->updateEr();
    }
    void setP(Vector3 t_p) {
        p3 = t_p;
        T = this->TFromP3();
        updateVwithP3T();
        updateEp();
    }
    void setV(Vector3 t_v) {
        v = t_v;
        updateP3TwithV();
        updateEp();
    }

    /**
     * Total cross section scattering with Target \p t.
     * Use current \p this->T as initial kinetic energy.
     * Return the statistical weight.
     */
    virtual Value totalCrossSection(const Target& t) const = 0;

    /**
     * Sample final state DM particle. The changes are made directly.
     * Return the statistical weight.
     */
    virtual Value scatter(const Target& t, RandomNumber<Value>& rn) = 0;

    /**
     * Generate a Event from the current state.
     * @param weight The statistical weight.
     * @param kmsec If \p true the quantities \p t and \p r will be in unit of second and
     * kilometer, respectively, otherwise remain natural unit with GeV = 1.
     */
    Event<Vector3, Value> toEvent(Value weight=1.0, bool kmsec=true) {
        if (kmsec) {
            return Event<Vector3, Value> {t/units::sec, T, r/units::km, p3, weight};
        }
        return Event<Vector3, Value> {t, T, r, p3, weight};
    }
};


/**
 * Base abstract Medium class.
 */
template<typename Vector3, typename Value>
class Medium
{
protected:
    std::vector<Target> targets;
public:
    Medium() {};
    Medium(const Medium&) = delete;
    Medium& operator=(const Medium&) = delete;
    Medium& operator=(const Medium&&) = delete;
    virtual ~Medium() {};

    virtual void setCache(const Particle<Vector3, Value>&) {}
    virtual const std::vector<Target>& getTargets() const { return targets; };
    /**
     * Propagate a particle.
     * Sample a free path and propagate the particle and determine whether
     * the particle is in the Medium or not.
     * Set \p Particle::in_medium to false if the simulation should stop.
     */
    virtual Value propagate(Particle<Vector3, Value>& particle,
                            RandomNumber<Value>& rn) = 0;
    /**
     * Sample a target before scattering.
     */
    virtual Target sampleTarget(const Particle<Vector3, Value>& particle,
                                RandomNumber<Value>& rn) const = 0;
};

/**
 * Basic interface for Particle scattering in the Medium.
 *
 * Update \p T, \p p3, \p v, \p ep of the particle \p p after a scattering.
 * It simply wraps \p Medium::sampleTarget and \p Particle::scatter.
 * Return the statistical weight.
 */
template<typename Vector3, typename Value>
Value scatter(Particle<Vector3, Value>& p,
              Medium<Vector3, Value>& m,
              RandomNumber<Value>& rn)
{
    if (p.in_medium) {
        Target target = m.sampleTarget(p, rn);
        return p.scatter(target, rn);
    }
    return 0.0;  // No scattering
}


/**
 * Basic interface for \p Particle propagating in the \p Medium.
 *
 * It simply wraps \p Medium::propagate. Return the statistical weight.
 */
template<typename Vector3, typename Value>
Value propagate(Particle<Vector3, Value>& p,
                       Medium<Vector3, Value>& m,
                       RandomNumber<Value>& rn)
{
    if (p.in_medium) {
        return m.propagate(p, rn);
    }
    return 0.0;  // No propagation
}

} // namespace darkprop
