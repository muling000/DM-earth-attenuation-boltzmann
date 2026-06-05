#include "Test.hpp"

TEST_CASE("halo velocity distribution") {
    REQUIRE(norm_maxwell() ==
            doctest::Approx(2.2005474759025309584e-9).epsilon(PREC).scale(REL));
    REQUIRE(norm_maxwell(constants::vesc) ==
            doctest::Approx(3.3270505493877547195e-8).epsilon(PREC).scale(REL));
    REQUIRE(norm_maxwell(0) == 0.0);

    REQUIRE(norm_esc() ==
            doctest::Approx(2.1859375381296047469e-9).epsilon(PREC).scale(REL));
    REQUIRE(norm_esc(1e100) ==
            doctest::Approx(norm_maxwell()).epsilon(PREC).scale(REL));

    REQUIRE(fv_halo_1d(0.0) == 0.0);
    REQUIRE(fv_halo_1d(constants::vesc) ==
            doctest::Approx(41.849179100293293953).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d(constants::v0) ==
            doctest::Approx(1138.8895021911499630).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d(7.9*units::km/units::sec) ==
            doctest::Approx(3.9868040069737837742).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d(7.9*units::km/units::sec,
                       600*units::km/units::sec, 240*units::km/units::sec) ==
               doctest::Approx(3.0690566193951478050).epsilon(PREC).scale(REL));

    REQUIRE(fv_halo_1d_earth(0.0) == 0.0);
    REQUIRE(fv_halo_1d_earth(constants::vesc+constants::vearth+PREC) == 0.0);
    REQUIRE(fv_halo_1d_earth(constants::vesc) ==
            doctest::Approx(256.04918983324484).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d_earth(constants::vesc - constants::vearth
                             + 1*units::km/units::sec) ==
            doctest::Approx(899.1758510342675).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d_earth(constants::vesc - constants::vearth
                             - 1*units::km/units::sec) ==
            doctest::Approx(897.9793090057444).epsilon(PREC).scale(REL));
    REQUIRE(fv_halo_1d_earth(7.9*units::km/units::sec, 240*units::km/units::sec,
                             600*units::km/units::sec, 250*units::km/units::sec) ==
            doctest::Approx(1.084802158967579).epsilon(PREC).scale(REL));
    CHECK_THROWS_AS(fv_halo_1d_earth(0.0, 0.0), std::invalid_argument);
}
