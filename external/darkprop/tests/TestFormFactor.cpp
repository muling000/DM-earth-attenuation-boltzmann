#include "Test.hpp"
#include "../darkprop/FormFactor.hpp"

#ifdef __cpp_lib_math_special_functions
Value cmath_core(Value x) {
    return 3 * std::sph_bessel(1, x) / x;
}

TEST_CASE("_ff_Helm_core") {
    auto xs = geomspace(1e-15, 1e3, 1000);

    for (Value x : xs) {
        CHECK(_ff_Helm_core(x) ==
              doctest::Approx(cmath_core(x)).epsilon(1e-2).scale(REL));
    }
}
#endif

TEST_CASE("ff_Helm zero input") {
    for (int a = 7; a < 200; ++a) {
        CHECK_EQ(ff_Helm(0.0, a), 1.0);
    }
}

TEST_CASE("ff_Helm exception A less than 7") {
    for (int a = 0; a < 7; ++a) {
        CHECK_THROWS_AS(ff_Helm(1.0, a), std::runtime_error);
    }
}
