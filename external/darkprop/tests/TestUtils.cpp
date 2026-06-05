#include "Test.hpp"

TEST_CASE("geomspace base10") {
    auto v = geomspace(1e-20, 1e20, 41);
    for (std::size_t i = 0; i != v.size(); ++i) {
        CHECK(v[i] == doctest::Approx(
                      1e-20*std::pow(10, i)).epsilon(PREC).scale(REL));
    }
}

TEST_CASE("linspace") {
    auto v = linspace(-1.0, 1.0, 11, true);
    for (std::size_t i = 0; i != v.size(); ++i) {
        CHECK(v[i] == doctest::Approx(-1.0 + i*0.2).epsilon(PREC).scale(ABS));
    }
}
