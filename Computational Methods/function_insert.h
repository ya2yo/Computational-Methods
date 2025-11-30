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
            throw invalid_argument("插值节点数量不匹配或为空");
        }
    }

    virtual ~FunctionInsert() {}


    virtual double interpolate(double val) = 0;


    void printNodes() const {
        cout << "插值节点：" << endl;
        for (int i = 0; i < n; i++) {
            cout << "(" << x[i] << ", " << y[i] << ") ";
        }
        cout << endl;
    }
};

// --- 拉格朗日插值 ---
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
            throw invalid_argument("赫米特插值要求节点、函数值和导数值数量一致");
        }
        computeHermiteCoeffs();
    }


    double interpolate(double val) override;
};

//三次样条函数插值 (采用边界条件 M0=Mn=0)
class CubicSplineInterpolation : public FunctionInsert {
private:
    vector<double> M; 

    void solveTridiagonalSystem();

public:
    CubicSplineInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : FunctionInsert(x_nodes, y_nodes) {
        M.resize(n, 0.0);
        if (n < 3) {
            throw invalid_argument("三次样条插值至少需要3个节点");
        }
        solveTridiagonalSystem(); 
    }
    double interpolate(double val) override;
    int findInterval(double val);
};