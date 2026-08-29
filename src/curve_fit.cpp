#include "curve_fit.h"

double Data::evaluate(double z) const {
    if (factor.empty()) throw runtime_error("Fit coefficients have not been computed");
    double result = 0.0;
    for (size_t i = 0; i < factor.size(); ++i) result += factor[i] * power(z, static_cast<int>(i));
    return result;
}

void Ordinary_Method::calculateFactor(int degree) {
    factor.clear();
    if (degree < 0 || static_cast<size_t>(degree) >= x.size()) {
        throw runtime_error("Polynomial degree must be between 0 and the number of data points minus one");
    }
    Matrix A(degree + 1, Vector(degree + 1, 0.0));
    Vector b(degree + 1, 0.0);
    int count = static_cast<int>(x.size());
    for (int i = 0; i <= degree; ++i) {
        for (int j = i; j <= degree; ++j) {
            double sum = 0.0;
            for (int k = 0; k < count; ++k) sum += power(x[k], i + j);
            A[i][j] = A[j][i] = sum;
        }
        for (int k = 0; k < count; ++k) b[i] += y[k] * power(x[k], i);
    }
    factor = solveLinearSystem(A, b);
}

void Ordinary_Method::calculate(int degree) { calculateFactor(degree); }

void Orthogonal_Method::calculateCoefficients(int degree) {
    factor.clear(); alpha.clear(); beta.clear();
    if (degree < 0 || static_cast<size_t>(degree) >= x.size()) {
        throw runtime_error("Polynomial degree must be between 0 and the number of data points minus one");
    }
    int count = static_cast<int>(x.size());
    factor.resize(degree + 1); alpha.resize(degree + 1, 0.0); beta.resize(degree + 1, 0.0);
    Vector current(count), previous(count, 0.0), previous_previous(count, 0.0);
    double sum_y_phi = 0.0, sum_phi_sq = static_cast<double>(count), sum_x_phi_sq = 0.0;
    for (int k = 0; k < count; ++k) { current[k] = 1.0; sum_y_phi += y[k]; sum_x_phi_sq += x[k]; }
    factor[0] = sum_y_phi / sum_phi_sq;
    alpha[0] = sum_x_phi_sq / sum_phi_sq;
    previous = current;
    double previous_sq = sum_phi_sq, previous_previous_sq = 0.0;
    for (int i = 1; i <= degree; ++i) {
        if (i >= 2) {
            if (previous_previous_sq < 1e-12) throw runtime_error("Orthogonal basis norm is too small");
            beta[i - 1] = previous_sq / previous_previous_sq;
        }
        sum_phi_sq = sum_y_phi = sum_x_phi_sq = 0.0;
        for (int k = 0; k < count; ++k) {
            double current_beta = i == 1 ? 0.0 : beta[i - 1];
            current[k] = (x[k] - alpha[i - 1]) * previous[k] - current_beta * previous_previous[k];
            sum_phi_sq += current[k] * current[k]; sum_y_phi += y[k] * current[k];
            sum_x_phi_sq += x[k] * current[k] * current[k];
        }
        if (sum_phi_sq < 1e-12) throw runtime_error("Orthogonal basis norm is too small");
        factor[i] = sum_y_phi / sum_phi_sq;
        if (i < degree) alpha[i] = sum_x_phi_sq / sum_phi_sq;
        previous_previous = previous; previous = current;
        previous_previous_sq = previous_sq; previous_sq = sum_phi_sq;
    }
}

void Orthogonal_Method::calculate(int degree) { calculateCoefficients(degree); }

double Orthogonal_Method::evaluate(double z) const {
    if (factor.empty()) throw runtime_error("Fit coefficients have not been computed");
    int degree = static_cast<int>(factor.size()) - 1;
    double previous_previous = 0.0, previous = 1.0, current = 0.0;
    double result = factor[0] * previous;
    for (int i = 1; i <= degree; ++i) {
        double current_beta = i == 1 ? 0.0 : beta[i - 1];
        current = (z - alpha[i - 1]) * previous - current_beta * previous_previous;
        result += factor[i] * current; previous_previous = previous; previous = current;
    }
    return result;
}
