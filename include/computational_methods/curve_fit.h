#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "direct_solvers.h"

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

inline double power(double x, int exponent) { return pow(x, exponent); }

inline Vector solveLinearSystem(const Matrix& A, const Vector& b) {
    if (A.empty() || b.empty() || A.size() != b.size()) {
        throw runtime_error("The matrix and vector must be non-empty and have matching dimensions");
    }
    try {
        CholeskyDecomposition solver(A, b);
        solver.solve();
        return solver.getSolution();
    } catch (const runtime_error&) {
        try {
            GaussianElimination solver(A, b);
            solver.solve();
            return solver.getSolution();
        } catch (const runtime_error& error) {
            throw runtime_error(string("Both Cholesky decomposition and Gaussian elimination failed: ") + error.what());
        }
    }
}

class Data {
protected:
    vector<double> x;
    vector<double> y;
    vector<double> factor;
public:
    Data(vector<double>& x_val, vector<double>& y_val) : x(x_val), y(y_val) {
        if (x.size() != y.size() || x.empty()) throw runtime_error("Data points must be non-empty and have matching x/y sizes");
    }
    virtual double evaluate(double z) const;
};

class Ordinary_Method : public Data {
private:
    void calculateFactor(int degree);
public:
    Ordinary_Method(vector<double>& x_val, vector<double>& y_val) : Data(x_val, y_val) {}
    void calculate(int degree);
    const vector<double>& getFactor() const { return factor; }
};

class Orthogonal_Method : public Data {
private:
    Vector alpha;
    Vector beta;
    void calculateCoefficients(int degree);
public:
    Orthogonal_Method(vector<double>& x_val, vector<double>& y_val) : Data(x_val, y_val) {}
    void calculate(int degree);
    double evaluate(double z) const override;
    const vector<double>& getFactor() const { return factor; }
};
