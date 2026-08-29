#include "numerical_integration.h"
#include <algorithm>

//¼ÆËãÌÝÐÎ¹«Ê½
double Romberg::trapezoid(int k) {
    double range = sup - inf;
    int n = 1 << k;
    double h = range / n;
    if (k == 0) {
        return range / 2.0 * (f(inf) + f(sup));
    }
    else {
        double sum_new = 0.0;
        int num_new_points = n / 2;
        for (int i = 1; i <= num_new_points; i++) {
            double x_i = inf + (2.0 * i - 1.0) * h;
            sum_new += f(x_i);
        }
        return 0.5 * R[k - 1][0] + h * sum_new;
    }
}

// ´òÓ¡ Romberg ±í¸ñ
void Romberg::printTable() const {
    cout << "\n--- Romberg »ý·Ö±í ---" << endl;
    for (int j = 0; j <= printn; j++) {
        cout << "R[" << j << "]: ";
        for (int i = 0; i <= j; i++) {
            cout << R[j][i] << " ";
        }
        cout << endl;
    }
    cout << "----------------------\n" << endl;
}
// Ö´ÐÐ Romberg »ý·Ö
double Romberg::solve(int max_k, double tolerance) {
    if (sup == inf) return 0.0;
    R.resize(max_k, vector<double>(max_k, 0.0));
    for (int k = 0; k < max_k; k++) {
        R[k][0] = trapezoid(k);
        for (int m = 1; m <= k; m++) {
            double factor = pow(4.0, m);
            R[k][m] = R[k][m - 1] + (R[k][m - 1] - R[k - 1][m - 1]) / (factor - 1.0);
        }
        if (k > 0) {
            double current_val = R[k][k];
            double prev_val = R[k - 1][k - 1];
            if (abs(current_val - prev_val) < tolerance) {
                printn = k;
                return current_val;
            }
        }
    }
    return R[max_k - 1][max_k - 1];
}