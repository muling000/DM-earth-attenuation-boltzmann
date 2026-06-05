/**
 * @file
 * The abstract MediumBall with spherical geometry.
 */

#pragma once

#include "Base.hpp"
#include "Vector3d.hpp"

namespace darkprop {

/**
 * Base class for spherical medium.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class MediumBall : public Medium<Vector3, Value>
{
protected:
    Value radius;  ///< the radius of the medium
    Value rfinal;  ///< a small distance that the particle travels after leaving the medium
public:
    explicit MediumBall(Value r=0.0, Value rf=0.0) : radius {r}, rfinal {rf} {};
    Value getRadius() const { return radius; }
    void setRadius(Value r) {
        if (r < 0) { throw std::runtime_error("radius less than 0"); }
        radius = r;
    }
    Value getRFinal() const { return rfinal; }
    void setRFinal(Value rf) {
        if (rf < 0) { throw std::runtime_error("rfinal less than 0"); }
        rfinal = rf;
    }
};

} // namespace darkprop
