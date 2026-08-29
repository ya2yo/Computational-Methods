#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace std;
// Direct solver base class
class DirectMethod {
protected:
    vector<vector<double>> m_A;  // Coefficient matrix
    vector<double> m_B;          // Right-hand side vector
    vector<double> ans;          // Solution vector
    int n;                       // Matrix dimension
public:
    DirectMethod(const vector<vector<double>>& A, const vector<double>& b)
        : m_A(A), m_B(b), n(static_cast<int>(b.size())) {
        ans.resize(n, 0);
    }
    virtual void solve() = 0;
    virtual ~DirectMethod() {}
    void print() const {
        cout << "The solution is:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = " << fixed << setprecision(6) << ans[i] << endl;
        }
    }
    vector<double> getSolution() const {
        return ans;
    }
};

// Gaussian elimination
class GaussianElimination : public DirectMethod {
public:
    GaussianElimination(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
    }
    void solve() override;
};

// LU decomposition
class LU_Decomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;  // Lower triangular matrix
    vector<vector<double>> m_U;  // Upper triangular matrix
    void doolittleDecompose();
public:
    LU_Decomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_U.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};

// Thomas algorithm for tridiagonal systems
class ChaseMethod : public DirectMethod {
private:
    vector<double> m_a;  // Subdiagonal
    vector<double> m_c;  // Superdiagonal
    vector<double> m_p;  // Decomposition auxiliary values
    vector<double> m_q;  // Decomposition auxiliary values
    bool decomposed;
    void decompose();
public:
    // Constructor for a tridiagonal matrix
    ChaseMethod(const vector<double>& a, const vector<double>& b,
        const vector<double>& c, const vector<double>& d);
    void solve() override;
};

// Cholesky decomposition
class CholeskyDecomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;   // Lower triangular matrix
    vector<vector<double>> m_LT;  // Transposed matrix
    void choleskyDecompose();
public:
    CholeskyDecomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_LT.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};