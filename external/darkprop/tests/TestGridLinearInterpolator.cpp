#include "Test.hpp"
#include <vector>
#include "../darkprop/numerical/GridLinearInterpolator.hpp"

using namespace numerical;

Value test_func_interp_1d(Value x)
{
    return std::exp(-x * x);
}

Value test_func_interp_2d(Value x, Value y)
{
    return x / (y * y + 1);
}

Value test_func_interp_3d(Value x, Value y, Value z)
{
    return x * x * x + y * y + z;
}

TEST_CASE("GridLinearInterpolator") {
    auto xs = linspace<Value>(-1, 1, 10);
    auto ys = linspace<Value>(-1, 1, 10);
    auto zs = linspace<Value>(-1, 1, 10);
    std::vector<Value> f_1d;
    std::vector<Value> f_2d;
    std::vector<Value> f_3d;
    for (auto x : xs) {
        f_1d.push_back(test_func_interp_1d(x));
        for (auto y : ys) {
            f_2d.push_back(test_func_interp_2d(x, y));
            for (auto z : zs) {
                f_3d.push_back(test_func_interp_3d(x, y, z));
            }
        }
    }

    Value fill_value = 1;
    auto interp_1d = GridLinearInterpolator(std::array{xs}, f_1d, fill_value);
    GridLinearInterpolator<2, std::vector<Value>> interp_2d;
    interp_2d.interpolate({xs, ys}, f_2d, fill_value);
    GridLinearInterpolator interp_3d {std::array{xs, ys, zs}, f_3d, fill_value};

    SUBCASE("PointsOnGrid") {
        for (auto x : xs) {
            CHECK(test_func_interp_1d(x) == doctest::Approx(
                  interp_1d.interp({x})).epsilon(PREC));
            for (auto y : ys) {
                CHECK(test_func_interp_2d(x, y) == doctest::Approx(
                      interp_2d.interp({x, y})).epsilon(PREC));
                for (auto z : zs) {
                    CHECK(test_func_interp_3d(x, y, z) == doctest::Approx(
                          interp_3d.interp({x, y, z})).epsilon(PREC));
                }
            }
        }
    }
    SUBCASE("PointsOffGrid") {
        auto testxs = linspace<Value>(-1, 1, 7);
        auto testys = linspace<Value>(-1, 1, 7);
        auto testzs = linspace<Value>(-1, 1, 7);
        for (auto x : testxs) {
            CHECK(test_func_interp_1d(x) == doctest::Approx(
                  interp_1d(x)).epsilon(0.1));
            for (auto y : testys) {
                CHECK(test_func_interp_2d(x, y) == doctest::Approx(
                      interp_2d(x, y)).epsilon(0.1));
                for (auto z : testzs) {
                    CHECK(test_func_interp_3d(x, y, z) == doctest::Approx(
                          interp_3d(x, y, z)).epsilon(0.1));
                }
            }
        }
    }
    SUBCASE("OutOfRange") {
        CHECK_EQ(interp_1d(-10.), fill_value);
        CHECK_EQ(interp_2d(-2., 0.), fill_value);
        CHECK_EQ(interp_3d(-2., 0., 0.), fill_value);
        CHECK_EQ(interp_3d(1., 1.5, -1.), fill_value);
        CHECK_EQ(interp_3d(0.5, 0.5, 3.0), fill_value);
    }
}
