#pragma once

#include <spdlog/spdlog.h>
#include "Const.hpp"
#include "Earth.hpp"
#include "SIDM.hpp"

namespace darkprop {

template<typename Vector3=Vector3d<double>, typename Value=double>
class HomoElectronEarth : public Earth<Vector3, Value>
{
private:
    static constexpr Value electron_density = 8e23 / (units::cm*units::cm*units::cm);
    const Target electron {"electron", -1, -1, constants::me};
    bool use_cache=false;
    Value cached_lambda = 0;

public:
    explicit HomoElectronEarth() : Earth<Vector3, Value>() {
        this->targets = std::vector<Target> { Target({"electron", -1, -1, constants::me}) };
    }
    ~HomoElectronEarth() {}
    Value meanFreePath(const Particle<Vector3, Value>& p) const {
        if (!(p.totalCrossSection(electron) > 0)) { // also checks nan
            spdlog::error("total cross section: {}", p.totalCrossSection(electron));
            throw std::runtime_error("total cross section is not grater than 0");
        }
        return 1.0 / (electron_density * p.totalCrossSection(electron));
    }
    void setCache(const Particle<Vector3, Value>& p) override {
        cached_lambda = meanFreePath(p);
        use_cache = true;
    }
    Value propagate(Particle<Vector3, Value>& p,
                    RandomNumber<Value>& rn) override {
        Value xi = rn.uniform_xi();
        Value freep {0.0};
        if (use_cache) {
            freep = -std::log(1 - xi) * cached_lambda;
        } else {
            freep = -std::log(1 - xi) * meanFreePath(p);
        }
        Vector3 rnow = p.r;
        p.r += p.ep * freep;
        if (p.r.norm() > this->radius) {
            p.in_medium = false;
            Value rnorm = rnow.norm();
            Value costheta = p.ep.dot(-p.er);
            freep = rnorm * costheta
                + std::sqrt(rnorm*rnorm*(costheta*costheta - 1) + this->radius*this->radius)
                + this->rfinal;
            p.r = rnow + p.ep * freep;
        }
        p.t += freep / p.v.norm();
        p.updateEr();
        return 1.0;
    }

    Target sampleTarget([[maybe_unused]] const Particle<Vector3, Value>& particle,
                        [[maybe_unused]] RandomNumber<Value>& rn) const override {
        return electron;
    }
};

} // namespace darkprop
