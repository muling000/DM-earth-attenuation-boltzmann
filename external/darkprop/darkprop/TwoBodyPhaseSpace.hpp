/**
 * @file
 * General sampling of two body phase space.
 */
#pragma once
#include <cmath>
#include "Base.hpp"
#include "Utils.hpp"


namespace darkprop::twobodyphase {


template<typename Value>
std::pair<Value, Value> E1_min_max(Value E0, Value p0_norm, Value m1, Value m2)
{
    Value s = E0 * E0 - p0_norm * p0_norm;
    Value term1 = (s + m1 * m1 - m2 * m2) * E0;
    Value term2 = p0_norm*std::sqrt((s-(m1+m2)*(m1+m2))*(s-(m1-m2)*(m1-m2)));

    Value E1Min = (term1 - term2) / (2.0 * s);
    Value E1Max = (term1 + term2) / (2.0 * s);

    return {E1Min, E1Max};
}


template<typename Value>
std::pair<Value, Value> sample_uniform(
    Value E0, Value p0_norm, Value m1, Value m2, RandomNumber<Value>& rn)
{
    const auto [E1min, E1max] = E1_min_max(E0, p0_norm, m1, m2);
    Value E1 = rn.uniform_ab(E1min, E1max);
    Value phi1 = rn.uniform_phi();
    return {E1, phi1};
}


template<typename Value>
Value costheta1(Value E0, Value p0_norm, Value m1, Value m2, Value E1)
{
    Value s = E0 * E0 - p0_norm * p0_norm;
    Value numerator = 2 * E0 * E1 - (s + m1 * m1 - m2 * m2);
    Value denominator = 2 * p0_norm * std::sqrt(E1 * E1 - m1 * m1);

    // TODO check the range [-1, 1]
    return numerator / denominator;
}


template<typename Vector3, typename Value>
Vector3 final_momentum(const Vector3& p0,
                       const Vector3& p1in,
                       Value p1out,
                       Value costheta1,
                       Value phi1)
{
    Vector3 ey = (p1in.cross(p0)).normalized();
    Vector3 ep0 = p0.normalized();
    Vector3 p1 = ep0 * p1out;
    rotation_axis_angle(p1, ey, std::acos(costheta1));
    rotation_axis_angle(p1, ep0, phi1);
    return p1;
}

} // namespace darkprop::twobodyphase
