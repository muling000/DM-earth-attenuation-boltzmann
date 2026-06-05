#include "Test.hpp"
#include <cmath>
#include <vector>
#include "../darkprop/numerical/Integrator.hpp"

using namespace numerical;

Value test_func_1d(double x)
{
    return x * x;
}

Value test_func_3d(const double *x, std::size_t dim)
{
    Value result = std::pow(2*M_PI, -1.5);
    for (std::size_t i = 0; i < dim; ++i) {
        result *= std::exp(-x[i]*x[i] / 2);
    }
    return result;
}

TEST_CASE("Integrator") {
    Integrator<Value (double)> integrator;
    CHECK(integrator(test_func_1d, 0, 1) == doctest::Approx(1./3).epsilon(1e-3));

    Integrator<Value (double)> integrator_log(true);
    CHECK(integrator_log(test_func_1d, 1e-15, 1) == doctest::Approx(1./3).epsilon(1e-3));

    SUBCASE("MCIntegrator") {
        MCIntegrator<Value (const double*, std::size_t)> mcintegrator;
        std::vector<double> xl {-10, -10, -10};
        std::vector<double> xu {10, 10, 10};
        CHECK(mcintegrator(test_func_3d, xl, xu) == doctest::Approx(1).epsilon(1e-3));
    }
}

