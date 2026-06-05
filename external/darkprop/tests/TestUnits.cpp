#include "Test.hpp"


TEST_CASE("unit consistency") {
    using namespace units;
    using namespace constants;
    CHECK(1.0 == doctest::Approx(299792458 * m / sec).epsilon(1e-7).scale(0));
    CHECK(1.0 == doctest::Approx(1.054571817e-34 * J * sec).epsilon(1e-7).scale(0));
    CHECK(1.0 == doctest::Approx(6.582119569e-22 * MeV * sec).epsilon(1e-7).scale(0));
    CHECK(me == doctest::Approx(9.1093837015e-31*kg).epsilon(1e-7).scale(0));
    CHECK(mp == doctest::Approx(1.67262192369e-27*kg).epsilon(1e-7).scale(0));
    CHECK(mn == doctest::Approx(1.00866491595*mu).epsilon(1e-7).scale(0));
    CHECK(mu == doctest::Approx(1.6605390666e-27*kg).epsilon(1e-7).scale(0));

    CHECK(F == doctest::Approx(C*C/J).epsilon(1e-9).scale(0));
    CHECK(8.987551e9*N == doctest::Approx(C*C/(4*M_PI*m*m)).epsilon(1e-7).scale(0));
    CHECK(1.0/137.035999084 == doctest::Approx(
            std::pow(1.602176634e-19*C, 2)/(4*M_PI)).epsilon(1e-7).scale(0));

    CHECK(year == doctest::Approx(365.242189*day).epsilon(1e-9).scale(0));
}
