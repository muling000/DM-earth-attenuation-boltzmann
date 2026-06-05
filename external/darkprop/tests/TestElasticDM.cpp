#include "Test.hpp"
#include "../darkprop/SIHelmDM.hpp"


TEST_CASE("SIHelmDM") {
    SIHelmDM<Vector3, double> dm0;
    SIHelmDM dm;
    dm.m = 1e-3;
    dm.sigma0 = 1e-28 * units::cm * units::cm;
    SUBCASE("initCrossSectionCDF") {
        dm.initCrossSectionCDF(1e-40, 1e2, oxygen, 100);
        // init twice
        dm.initCrossSectionCDF(1e-50, 1e5, oxygen, 1000);
    }
    SUBCASE("inverseCDF") {
        dm.T = 1e-3;
        dm.initCrossSectionCDF(1e-50, 1e5, oxygen, 10000);
        dm.initCrossSectionCDF(1e-50, 1e5, potassium, 10000);
        CHECK_EQ(1e-50, dm.inverseCDF(0, oxygen));
        dm.inverseCDF(1e-15, oxygen);
        dm.inverseCDF(1e-15, potassium);
        dm.inverseCDF(0.5, oxygen);
        dm.inverseCDF(0.5, potassium);
        CHECK(dm.maximumRecoilT(oxygen) == doctest::Approx(
                    dm.inverseCDF(1, oxygen)).epsilon(PREC).scale(REL));
    }

}
