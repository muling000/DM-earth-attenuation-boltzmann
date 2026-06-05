/**
 * @file
 * The abstract Earth with radius \p constants::rEarth.
 */

#pragma once

#include "MediumBall.hpp"
#include "Const.hpp"

namespace darkprop {

/**
 * Abstract Earth class.
 */
template<typename Vector3=Vector3d<double>, typename Value=double>
class Earth : public MediumBall<Vector3, Value>
{
public:
    explicit Earth(Value r=constants::rEarth, Value rf=1.0*units::cm)
        : MediumBall<Vector3, Value>(r, rf) {};
};

} // namespace darkprop
