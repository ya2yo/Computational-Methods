#include "Direct.h"
#include "function_insert.h"
#include "linear_equations_iterative.h"
#include "numerical_integration.h"

#include <cassert>
#include <cmath>
#include <vector>

namespace {
constexpr double epsilon = 1e-8;

void expect_near(double actual, double expected) {
    assert(std::abs(actual - expected) < epsilon);
}

void test_gaussian_elimination() {
    const std::vector<std::vector<double>> matrix{{2.0, 1.0}, {1.0, 3.0}};
    const std::vector<double> rhs{5.0, 6.0};
    GaussianElimination solver(matrix, rhs);
    solver.solve();
    const auto result = solver.getSolution();
    expect_near(result[0], 1.8);
    expect_near(result[1], 1.4);
}

void test_interpolation() {
    const std::vector<double> x{0.0, 1.0, 2.0};
    const std::vector<double> y{0.0, 1.0, 4.0};
    LagrangeInterpolation interpolation(x, y);
    expect_near(interpolation.interpolate(1.5), 2.25);
}

void test_romberg() {
    Romberg integration(1.0, 0.0, [](double x) { return 4.0 / (1.0 + x * x); });
    expect_near(integration.solve(8, 1e-8), 3.14159265358979323846);
}
}

int main() {
    test_gaussian_elimination();
    test_interpolation();
    test_romberg();
    return 0;
}
