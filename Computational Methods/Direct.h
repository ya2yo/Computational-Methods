#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace std;
// 直接法基类
class DirectMethod {
protected:
    vector<vector<double>> m_A;  // 系数矩阵
    vector<double> m_B;          // 右端向量
    vector<double> ans;          // 解向量
    int n;                       // 矩阵维度
public:
    DirectMethod(const vector<vector<double>>& A, const vector<double>& b)
        : m_A(A), m_B(b), n(b.size()) {
        ans.resize(n, 0);
    }
    virtual void solve() = 0;
    virtual ~DirectMethod() {}
    void print() const {
        cout << "方程组的解为：" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = " << fixed << setprecision(6) << ans[i] << endl;
        }
    }
    vector<double> getSolution() const {
        return ans;
    }
};

// 高斯消去法
class GaussianElimination : public DirectMethod {
public:
    GaussianElimination(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
    }
    void solve() override;
};

// LU分解法
class LU_Decomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;  // 下三角矩阵
    vector<vector<double>> m_U;  // 上三角矩阵
    void doolittleDecompose();
public:
    LU_Decomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_U.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};

// 追赶法（三对角方程组）
class ChaseMethod : public DirectMethod {
private:
    vector<double> m_a;  // 下对角线
    vector<double> m_c;  // 上对角线
    vector<double> m_p;  // 分解中间变量
    vector<double> m_q;  // 分解中间变量
    bool decomposed;
    void decompose();
public:
    // 特殊构造函数用于三对角矩阵
    ChaseMethod(const vector<double>& a, const vector<double>& b,
        const vector<double>& c, const vector<double>& d);
    void solve() override;
};

// Cholesky分解法
class CholeskyDecomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;   // 下三角矩阵
    vector<vector<double>> m_LT;  // 转置矩阵
    void choleskyDecompose();
public:
    CholeskyDecomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_LT.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};