// #include <Eigen/Geometry>
#include "doctest.h"
#include "../darkprop/core.hpp"

using namespace darkprop;

// typedef Eigen::Vector3d Vector3;
typedef Vector3d<double> Vector3;
typedef double Value;

inline const Value PREC = std::numeric_limits<Value>::epsilon() * 100;
// const double PREC = 1e-15;
// inline const double PREC = Eigen::NumTraits<double>::dummy_precision();
inline constexpr Value ABS = 1.0;
inline constexpr Value REL = 0.0;
inline const Value sqrt2 = std::sqrt(2.0);

const Target oxygen {"O", 8, 16, 16 * constants::mu};
const Target potassium {"K", 19, 36, 36 * constants::mu};
