#include "interpolation.h"
#include "numerical_integration.h"
#include "iterative_solvers.h"
#include "direct_solvers.h"
#include "root_finding.h"
#include "ordinary_differential_equations.h"
#include "curve_fit.h"
double f2(double x) {
    return 0.5 * pow(10.0 - pow(x, 3), 0.5);
}

int main() {
    vector<vector<double>> A = {
        {2,0,0,0},
        {0.5,2,0.5,0},
        {0,0.5,2,0.5},
        {0,0,0,2}
    };
    vector<double> b = { -0.6,3,-6,6.6};
    LU_Decomposition ans(A, b);
    ans.solve();
    cout << "----------------------" << endl;
    cout << "LU decomposition" << endl;
    ans.print();
    cout << "----------------------" << endl;
    cout << "Gauss-Seidel method" << endl;
    Gauss_Seidel solution(A, b, 1e-5);
    solution.solve();
    solution.print();
    {
        auto f = [](double x) {
            return x * x * x - x * x - x - 1;
            };
        auto g = [](double x) {
            return 1 + 1.0 / x + 1.0 / x / x;
            };
        auto df = [](double x) {
            return 3.0 * x * x - 2.0 * x - 1;
            };
        cout << "----------------" << endl;
        cout << "Fixed-point iteration" << endl;
        Stationary_Points_Iterate a(f, g);
        a.solve(1);
        cout << "Newton method" << endl;
        Newton c(f, df);
        c.solve(2);
        cout << "Steffensen acceleration" << endl;
        Steffensen b(f, g);
        b.solve(1);
    }
    {
        cout << "-------------------" << endl;
        cout << "Curve fitting test" << endl;
        vector<double> x{ 1,3,4,6,7 };
        vector<double> y{ -2.1,-0.9,-0.6,0.6,0.9 };
        Ordinary_Method solution(x, y);
        solution.calculate(1);
        auto ans = solution.getFactor();
        cout << ans[0] << "+" << ans[1] << "x" << endl;
    }
    {
        cout << "-------------------" << endl;
        cout << "Romberg integration test" << endl;
        auto f = [](double x) {
            return 4.0 / (1 + x * x);
            };
        Romberg solution(1, 0, f);
        solution.solve();
        solution.printTable();
    }
    {
        cout << "----------------------" << endl;
        cout << "ODE numerical solution test" << endl;
        auto dy = [](double x, double y)->double {
            return x - y + 1;
            };
        Runge_Kutta4 solution(0.0,1.0, dy, 0.1, 0.5);
        solution.solve();
    }
}
