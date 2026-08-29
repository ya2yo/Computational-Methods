#include "function_insert.h"

// Compute the Lagrange basis polynomial L_k(x).
double lagrangeBasis(int k, double val, const vector<double>& x_nodes) {
    double result = 1.0;
    for (size_t i = 0; i < x_nodes.size(); ++i) {
        if (i != static_cast<size_t>(k)) {
            result *= (val - x_nodes[i]) / (x_nodes[k] - x_nodes[i]);
        }
    }
    return result;
}

double LagrangeInterpolation::interpolate(double val) {
    double polynomial = 0.0;
    for (int k = 0; k < n; ++k) polynomial += y[k] * lagrangeBasis(k, val, x);
    return polynomial;
}

void NewtonInterpolation::computeDividedDiffs() {
    divided_diffs = y;
    for (int k = 1; k < n; ++k) {
        for (int i = n - 1; i >= k; --i) {
            divided_diffs[i] = (divided_diffs[i] - divided_diffs[i - 1]) / (x[i] - x[i - k]);
        }
    }
}

double NewtonInterpolation::interpolate(double val) {
    double polynomial = divided_diffs[0];
    double product = 1.0;
    for (int k = 1; k < n; ++k) {
        product *= val - x[k - 1];
        polynomial += divided_diffs[k] * product;
    }
    return polynomial;
}

void HermiteInterpolation::computeHermiteCoeffs() {
    vector<double> z(2 * n), f(2 * n);
    for (int i = 0; i < n; ++i) {
        z[2 * i] = z[2 * i + 1] = x[i];
        f[2 * i] = f[2 * i + 1] = y[i];
    }
    H_coeffs = f;
    for (int i = 0; i < n; ++i) {
        if (i == 0) H_coeffs[1] = y_prime[i];
        else {
            H_coeffs[2 * i] = (H_coeffs[2 * i] - H_coeffs[2 * i - 1]) / (z[2 * i] - z[2 * i - 1]);
            H_coeffs[2 * i + 1] = y_prime[i];
        }
    }
    for (int k = 2; k < 2 * n; ++k) {
        for (int i = 2 * n - 1; i >= k; --i) {
            H_coeffs[i] = (H_coeffs[i] - H_coeffs[i - 1]) / (z[i] - z[i - k]);
        }
    }
}

double HermiteInterpolation::interpolate(double val) {
    vector<double> z(2 * n);
    for (int i = 0; i < n; ++i) z[2 * i] = z[2 * i + 1] = x[i];
    double value = H_coeffs[0], product = 1.0;
    for (int k = 1; k < 2 * n; ++k) {
        product *= val - z[k - 1];
        value += H_coeffs[k] * product;
    }
    return value;
}

void CubicSplineInterpolation::solveTridiagonalSystem() {
    if (n < 3) return;
    int dim = n - 2;
    vector<double> h(n - 1), a(dim), b(dim), c(dim), d(dim);
    for (int i = 0; i < n - 1; ++i) h[i] = x[i + 1] - x[i];
    for (int i = 0; i < dim; ++i) {
        int k = i + 1;
        if (k > 0) a[i] = h[k - 1];
        b[i] = 2.0 * (h[k - 1] + h[k]);
        if (k < n - 2) c[i] = h[k];
        d[i] = 6.0 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1]);
    }
    vector<double> p(dim), q(dim);
    p[0] = c[0] / b[0];
    q[0] = d[0] / b[0];
    for (int i = 1; i < dim; ++i) {
        double denominator = b[i] - a[i] * p[i - 1];
        if (abs(denominator) < 1e-10) throw runtime_error("The spline matrix is singular");
        if (i < dim - 1) p[i] = c[i] / denominator;
        q[i] = (d[i] - a[i] * q[i - 1]) / denominator;
    }
    vector<double> temp(dim);
    temp[dim - 1] = q[dim - 1];
    for (int i = dim - 2; i >= 0; --i) temp[i] = q[i] - p[i] * temp[i + 1];
    for (int i = 0; i < dim; ++i) M[i + 1] = temp[i];
}

int CubicSplineInterpolation::findInterval(double val) {
    if (val < x[0] || val > x[n - 1]) throw out_of_range("Interpolation point is outside the node range");
    for (int i = 0; i < n - 1; ++i) if (val >= x[i] && val <= x[i + 1]) return i;
    return n - 2;
}

double CubicSplineInterpolation::interpolate(double val) {
    if (n < 3) {
        if (n == 2) return y[0] + (y[1] - y[0]) * (val - x[0]) / (x[1] - x[0]);
        if (n == 1) return y[0];
        throw invalid_argument("At least one interpolation node is required");
    }
    int i = findInterval(val), next = i + 1;
    double h_i = x[next] - x[i];
    double u = (x[next] - val) / h_i, v = (val - x[i]) / h_i;
    return M[i] * u * (u * u - 1.0) * (h_i * h_i / 6.0) +
           M[next] * v * (v * v - 1.0) * (h_i * h_i / 6.0) + y[i] * u + y[next] * v;
}
