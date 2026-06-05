/**
 * @file
 * Tools to analyze the trajectories.
 */

#pragma once

#include <filesystem>
#include "Base.hpp"

namespace darkprop {

template<typename Vector3, typename Value>
inline void _push_back_event(const Event<Vector3, Value>& p1,
                             const Event<Vector3, Value>& p2,
                             Value l,
                             std::vector<Event<Vector3, Value>>& events)
{
    if (l >= 0.0 && l <= 1.0) {
        events.push_back(Event<Vector3, Value>(
                    (p1.t + l*(p2.t - p1.t))/units::sec,
                     p1.T,
                     (p1.r + l*(p2.r - p1.r))/units::km,
                     p1.p3,
                     p2.weight));
    }
}

template<typename Vector3, typename Value>
void sphere_cross(const Vector3& r0, Value R,
                  const Event<Vector3, Value>& p1,
                  const Event<Vector3, Value>& p2,
                  std::vector<Event<Vector3, Value>>& events)
{
    bool p1_in_sphere = (p1.r - r0).norm() < R;
    bool p2_in_sphere = (p2.r - r0).norm() < R;
    if (p1_in_sphere && p2_in_sphere) {
        return;
    }
    Vector3 r21 = p2.r - p1.r;
    Vector3 r10 = p1.r - r0;
    Value a = r21.dot(r21);
    Value b = 2. * r21.dot(r10);
    Value c = r10.dot(r10) - R*R;
    Value Delta = b*b - 4.*a*c;
    if (!p1_in_sphere && !p2_in_sphere) {
        if (Delta > 0) {
            Value l1 = (-b - std::sqrt(Delta)) / (2. * a);
            Value l2 = (-b + std::sqrt(Delta)) / (2. * a);
            _push_back_event(p1, p2, l1, events);
            _push_back_event(p1, p2, l2, events);
        } else if (Delta == 0) {
            Value l = -b / (2.*a);
            _push_back_event(p1, p2, l, events);
        }
    } else {
        Value l;
        if (p1_in_sphere) {
            l = (-b + std::sqrt(Delta)) / (2.*a);
        } else {
            l = (-b - std::sqrt(Delta)) / (2.*a);
        }
        _push_back_event(p1, p2, l, events);
    }
}

} // namespace darkprop
