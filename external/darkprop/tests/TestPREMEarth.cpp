#include "Test.hpp"

TEST_CASE("PREMEarth") {
    PREMEarth<Vector3, Value> prem;
    CHECK(prem.densityIntegral(0, 0, 0.0, 0) == 0);
    // TODO check for L != 0
    for (int layer = 1; layer < 10; ++layer) {
        Value r0 = prem.getLayerRadius(layer - 1);
        CHECK(prem.densityIntegral(r0, 0, 0, layer)
              == doctest::Approx(0).epsilon(PREC));
        CHECK(prem.densityIntegral(r0, 0.5, 0, layer)
              == doctest::Approx(0).epsilon(PREC));
        CHECK(prem.densityIntegral(r0, -1, 0, layer)
              == doctest::Approx(0).epsilon(PREC));
        CHECK(prem.densityIntegral(r0, 1, 0, layer)
              == doctest::Approx(0).epsilon(PREC));
    }
}
