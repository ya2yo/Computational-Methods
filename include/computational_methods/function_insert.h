#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace std;

class FunctionInsert {
protected:
    vector<double> x; 
    vector<double> y; 
    int n;            
public:

    FunctionInsert(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : x(x_nodes), y(y_nodes), n(x_nodes.size()) {
        if (x.size() != y.size() || x.empty()) {
            throw invalid_argument("Interpolation nodes and values must be non-empty and have matching sizes");
        }
    }

    virtual ~FunctionInsert() {}


    virtual double interpolate(double val) = 0;


    void printNodes() const {
        cout << "Interpolation nodes: " << endl;
        for (int i = 0; i < n; i++) {
            cout << "(" << x[i] << ", " << y[i] << ") ";
        }
        cout << endl;
    }
};

// --- Lagrange interpolation ---
class LagrangeInterpolation : public FunctionInsert {
public:
    LagrangeInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : FunctionInsert(x_nodes, y_nodes) {
    }


    double interpolate(double val) override;
};


class NewtonInterpolation : public FunctionInsert {
private:
    vector<double> divided_diffs;
    void computeDividedDiffs();

public:
    NewtonInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : FunctionInsert(x_nodes, y_nodes) {
        computeDividedDiffs(); 
    }


    double interpolate(double val) override;
};



class HermiteInterpolation : public FunctionInsert {
private:
    vector<double> y_prime;
    vector<double> H_coeffs; 

    void computeHermiteCoeffs(); 

public:
    HermiteInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes, const vector<double>& y_prime_nodes)
        : FunctionInsert(x_nodes, y_nodes), y_prime(y_prime_nodes) {
        if (x.size() != y_prime.size()) {
            throw invalid_argument("Hermite interpolation requires matching node, value, and derivative counts");
        }
        computeHermiteCoeffs();
    }


    double interpolate(double val) override;
};

// Cubic spline interpolation with natural boundary conditions
class CubicSplineInterpolation : public FunctionInsert {
private:
    vector<double> M; 

    void solveTridiagonalSystem();

public:
    CubicSplineInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : FunctionInsert(x_nodes, y_nodes) {
        M.resize(n, 0.0);
        if (n < 3) {
            throw invalid_argument("Cubic spline interpolation requires at least three nodes");
        }
        solveTridiagonalSystem(); 
    }
    double interpolate(double val) override;
    int findInterval(double val);
};