#include "function_insert.h"
#include "numerical_integration.h"
#include "linear_equations_iterative.h"
#include "Direct.h"
#include "Unary_Equation.h"
#include "Ordinary_Differential_Equations.h"
#include "CurveFit.h"
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
    cout << "LU·Ö½â·¨" << endl;
    ans.print();
    cout << "----------------------" << endl;
    cout << "¸ßË¹ÈüµÂ·¨" << endl;
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
        cout << "²»¶¯µã·¨" << endl;
        Stationary_Points_Iterate a(f, g);
        a.solve(1);
        cout << "Å£¶Ù·¨" << endl;
        Newton c(f, df);
        c.solve(2);
        cout << "Steffensen¼ÓËÙ·¨" << endl;
        Steffensen b(f, g);
        b.solve(1);
    }
    {
        cout << "-------------------" << endl;
        cout << "ÇúÏßÄâºÏ²âÊÔ" << endl;
        vector<double> x{ 1,3,4,6,7 };
        vector<double> y{ -2.1,-0.9,-0.6,0.6,0.9 };
        Ordinary_Method solution(x, y);
        solution.calculate(1);
        auto ans = solution.getFactor();
        cout << ans[0] << "+" << ans[1] << "x" << endl;
    }
    {
        cout << "-------------------" << endl;
        cout << "ÊýÖµ»ý·ÖRomberg²âÊÔ" << endl;
        auto f = [](double x) {
            return 4.0 / (1 + x * x);
            };
        Romberg solution(1, 0, f);
        solution.solve();
        solution.printTable();
    }
    {
        cout << "----------------------" << endl;
        cout << "³£Î¢·Ö·½³ÌÊýÖµ½â·¨²âÊÔ" << endl;
        auto dy = [](double x, double y)->double {
            return x - y + 1;
            };
        Runge_Kutta4 solution(0.0,1.0, dy, 0.1, 0.5);
        solution.solve();
    }
}
