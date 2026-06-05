#include "Test.hpp"

TEST_CASE("inject") {
    HomoEarth<Vector3, Value> earth;
    SIDM<Vector3, Value> dm;
    CHECK_FALSE(dm.in_medium);

    SUBCASE("ordinary") {
        dm.setP({0, 0, -1});
        dm.setR({0, 0, constants::rEarth * 1.1});
        inject(dm);
        CHECK(dm.in_medium);
        dm.in_medium = false;

        dm.setR({0, 0, constants::rEarth});
        inject(dm);
        CHECK(dm.in_medium);
        dm.in_medium = false;

        dm.setR({0, 0, constants::rEarth * 0.9});
        inject(dm);
        CHECK(dm.in_medium);
        dm.in_medium = false;

        dm.setR({constants::rEarth / 2, 0, constants::rEarth * 2});
        inject(dm);
        CHECK(dm.in_medium);
    }

    SUBCASE("exception") {
        dm.setP({0, 0, 1});
        dm.setR({0, 0, constants::rEarth * 1.1});
        CHECK_THROWS_AS(inject(dm), std::runtime_error);
        dm.setP({-1, 1, 0});
        CHECK_THROWS_AS(inject(dm), std::runtime_error);
        dm.setP({0, 0, -10});
        dm.setR({constants::rEarth, 0, constants::rEarth * 1.1});
        CHECK_THROWS_AS(inject(dm), std::runtime_error);
        dm.setP({1, 0, 1});
        dm.setR({-constants::rEarth * sqrt2, 0, 0});
        CHECK_THROWS_AS(inject(dm), std::runtime_error);
    }

    SUBCASE("special") {
        dm.setP({0, 0, -10});
        dm.setR({constants::rEarth-1e-5*units::mm, 0, constants::rEarth*1.1});
        inject(dm);
        CHECK(dm.in_medium);
    }
}
