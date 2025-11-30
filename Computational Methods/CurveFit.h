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
        throw runtime_error("求解线性方程组时矩阵或向量大小不匹配或为空");
    }

    try {
        CholeskyDecomposition solver(A, b);
        solver.solve();
        return solver.getSolution();
    }
    catch (const runtime_error& e) {
        cerr << "警告: Cholesky 分解失败 (" << e.what() << ")，尝试使用高斯消元法求解..." << endl;
        try {
            GaussianElimination solver(A, b);
            solver.solve();
            return solver.getSolution();
        }
        catch (const runtime_error& ge) {
            throw runtime_error("高斯消元法也失败了: " + string(ge.what()));
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
            throw runtime_error("数据点的 x 和 y 个数不相等或数据为空!");
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