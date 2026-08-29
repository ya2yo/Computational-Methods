#include "curve_fit.h"
#include "direct_solvers.h"
#include "interpolation.h"
#include "iterative_solvers.h"
#include "numerical_integration.h"
#include "ordinary_differential_equations.h"
#include "root_finding.h"

#include <cassert>
#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double epsilon = 1e-7;

void expect_near(double actual, double expected, double tolerance = epsilon) {
    assert(std::abs(actual - expected) < tolerance);
}

void expect_vector_near(const std::vector<double>& actual, const std::vector<double>& expected,
                        double tolerance = epsilon) {
    assert(actual.size() == expected.size());
    for (size_t i = 0; i < actual.size(); ++i) expect_near(actual[i], expected[i], tolerance);
}

template <typename Exception, typename Callable>
void expect_throw(Callable&& callable) {
    bool thrown = false;
    try {
        callable();
    } catch (const Exception&) {
        thrown = true;
    }
    assert(thrown);
}

void test_direct_solvers() {
    const std::vector<std::vector<double>> matrix{{2.0, 1.0}, {1.0, 3.0}};
    const std::vector<double> rhs{5.0, 6.0};
    const std::vector<double> expected{1.8, 1.4};

    GaussianElimination gaussian(matrix, rhs);
    gaussian.solve();
    expect_vector_near(gaussian.getSolution(), expected);

    LU_Decomposition lu(matrix, rhs);
    lu.solve();
    expect_vector_near(lu.getSolution(), expected);

    CholeskyDecomposition cholesky(matrix, rhs);
    cholesky.solve();
    expect_vector_near(cholesky.getSolution(), expected);

    ChaseMethod thomas({0.0, 1.0}, {2.0, 2.0}, {1.0, 0.0}, {5.0, 6.0});
    thomas.solve();
    expect_vector_near(thomas.getSolution(), {4.0 / 3.0, 7.0 / 3.0});

    expect_throw<std::runtime_error>([] {
        GaussianElimination solver({{1.0, 2.0}, {2.0, 4.0}}, {3.0, 6.0});
        solver.solve();
    });
}

void test_iterative_solvers() {
    const std::vector<std::vector<double>> matrix{{4.0, 1.0}, {1.0, 3.0}};
    const std::vector<double> rhs{5.0, 6.0};
    const std::vector<double> expected{0.818181818, 1.727272727};

    Jacobi jacobi(matrix, rhs, 1e-8);
    jacobi.solve();
    expect_vector_near(jacobi.getSolution(), expected, 1e-5);

    Gauss_Seidel gauss_seidel(matrix, rhs, 1e-8);
    gauss_seidel.solve();
    expect_vector_near(gauss_seidel.getSolution(), expected, 1e-5);

    SOR sor(matrix, rhs, 1e-8, 1.1);
    sor.solve();
    expect_vector_near(sor.getSolution(), expected, 1e-5);
}

void test_interpolation() {
    const std::vector<double> x{0.0, 1.0, 2.0};
    const std::vector<double> y{0.0, 1.0, 4.0};
    LagrangeInterpolation lagrange(x, y);
    NewtonInterpolation newton(x, y);
    expect_near(lagrange.interpolate(1.5), 2.25);
    expect_near(newton.interpolate(1.5), 2.25);

    HermiteInterpolation hermite({0.0, 1.0}, {0.0, 1.0}, {1.0, 1.0});
    expect_near(hermite.interpolate(0.5), 0.0);

    CubicSplineInterpolation spline(x, y);
    expect_near(spline.interpolate(1.5), 2.3125);
    expect_throw<std::out_of_range>([&] { spline.interpolate(3.0); });
    expect_throw<std::invalid_argument>([] {
        LagrangeInterpolation invalid({0.0}, {0.0, 1.0});
        (void)invalid;
    });
}

void test_curve_fitting() {
    std::vector<double> x{0.0, 1.0, 2.0};
    std::vector<double> y{1.0, 3.0, 5.0};
    Ordinary_Method ordinary(x, y);
    ordinary.calculate(1);
    expect_near(ordinary.evaluate(3.0), 7.0, 1e-5);

    Orthogonal_Method orthogonal(x, y);
    orthogonal.calculate(1);
    expect_near(orthogonal.evaluate(3.0), 7.0, 1e-5);
    expect_throw<std::runtime_error>([&] { ordinary.calculate(3); });
}

double root_function(double x) { return x * x - 2.0; }
double root_derivative(double x) { return 2.0 * x; }
double fixed_point(double x) { return 0.5 * (x + 2.0 / x); }

double ode_function(double, double y) { return y; }

void test_root_finding_and_ode() {
    Newton newton(root_function, root_derivative);
    std::istringstream input("1\\n");
    auto* original_buffer = std::cin.rdbuf(input.rdbuf());
    const double newton_root = newton.solve(1.5);
    std::cin.rdbuf(original_buffer);
    expect_near(newton_root, std::sqrt(2.0), 1e-5);

    Stationary_Points_Iterate fixed(root_function, fixed_point);
    expect_near(fixed.solve(1.5), std::sqrt(2.0), 1e-5);

    Steffensen steffensen(root_function, fixed_point);
    expect_near(steffensen.solve(1.5), std::sqrt(2.0), 1e-5);

    Runge_Kutta4 ode(0.0, 1.0, ode_function, 0.1, 0.2);
    ode.solve();
    const auto values = ode.gety();
    expect_near(values.back(), std::exp(0.2), 1e-4);
}

void test_romberg() {
    Romberg integration(1.0, 0.0, [](double x) { return 4.0 / (1.0 + x * x); });
    expect_near(integration.solve(8, 1e-8), 3.14159265358979323846, 1e-6);
    Romberg reversed(0.0, 1.0, [](double x) { return x; });
    expect_near(reversed.solve(), 0.5);
    Romberg constant(1.0, 1.0, [](double) { return 42.0; });
    expect_near(constant.solve(), 0.0);
}
}

int main() {
    test_direct_solvers();
    test_iterative_solvers();
    test_interpolation();
    test_curve_fitting();
    test_root_finding_and_ode();
    test_romberg();
    return 0;
}
