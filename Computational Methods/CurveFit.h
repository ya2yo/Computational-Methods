#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "Direct.h"

using namespace std;


using Matrix = vector<vector<double>>;
using Vector = vector<double>;

inline double power(double x, int m) {
    return pow(x, m);
}
inline Vector solveLinearSystem(const Matrix& A, const Vector& b) {
    if (A.empty() || b.empty() || A.size() != b.size()) {
        throw runtime_error("Çó½âÏßÐÔ·½³Ì×éÊ±¾ØÕó»òÏòÁ¿´óÐ¡²»Æ¥Åä»òÎª¿Õ");
    }

    try {
        CholeskyDecomposition solver(A, b);
        solver.solve();
        return solver.getSolution();
    }
    catch (const runtime_error& e) {
        cerr << "¾¯¸æ: Cholesky ·Ö½âÊ§°Ü (" << e.what() << ")£¬³¢ÊÔÊ¹ÓÃ¸ßË¹ÏûÔª·¨Çó½â..." << endl;
        try {
            GaussianElimination solver(A, b);
            solver.solve();
            return solver.getSolution();
        }
        catch (const runtime_error& ge) {
            throw runtime_error("¸ßË¹ÏûÔª·¨Ò²Ê§°ÜÁË: " + string(ge.what()));
        }
    }
}
class Data {
protected:
    vector<double> x;
    vector<double> y;
    vector<double> factor; 
public:
    Data(vector<double>& x_val, vector<double>& y_val) :x(x_val), y(y_val) {
        if (x.size() != y.size() || x.empty()) {
            throw runtime_error("Êý¾ÝµãµÄ x ºÍ y ¸öÊý²»ÏàµÈ»òÊý¾ÝÎª¿Õ!");
        }
    }
    virtual double evaluate(double z) const;
};


class Ordinary_Method :public Data {
private:
    void calculateFactor(int n);
public:
    Ordinary_Method(vector<double>& x_val, vector<double>& y_val) :Data(x_val, y_val) {}
    void calculate(int n);
    const vector<double>& getFactor() const { return factor; }
};


class Orthogonal_Method : public Data {
private:
    Vector alpha;
    Vector beta;
    void calculateCoefficients(int n);

public:
    Orthogonal_Method(vector<double>& x_val, vector<double>& y_val) :Data(x_val, y_val) {}

    void calculate(int n);

    double evaluate(double z) const override;

    const vector<double>& getFactor() const { return factor; }
};