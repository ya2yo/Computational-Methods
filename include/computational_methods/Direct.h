#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace std;
// Ö±½Ó·¨»ùÀà
class DirectMethod {
protected:
    vector<vector<double>> m_A;  // ÏµÊý¾ØÕó
    vector<double> m_B;          // ÓÒ¶ËÏòÁ¿
    vector<double> ans;          // ½âÏòÁ¿
    int n;                       // ¾ØÕóÎ¬¶È
public:
    DirectMethod(const vector<vector<double>>& A, const vector<double>& b)
        : m_A(A), m_B(b), n(b.size()) {
        ans.resize(n, 0);
    }
    virtual void solve() = 0;
    virtual ~DirectMethod() {}
    void print() const {
        cout << "·½³Ì×éµÄ½âÎª£º" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = " << fixed << setprecision(6) << ans[i] << endl;
        }
    }
    vector<double> getSolution() const {
        return ans;
    }
};

// ¸ßË¹ÏûÈ¥·¨
class GaussianElimination : public DirectMethod {
public:
    GaussianElimination(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
    }
    void solve() override;
};

// LU·Ö½â·¨
class LU_Decomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;  // ÏÂÈý½Ç¾ØÕó
    vector<vector<double>> m_U;  // ÉÏÈý½Ç¾ØÕó
    void doolittleDecompose();
public:
    LU_Decomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_U.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};

// ×·¸Ï·¨£¨Èý¶Ô½Ç·½³Ì×é£©
class ChaseMethod : public DirectMethod {
private:
    vector<double> m_a;  // ÏÂ¶Ô½ÇÏß
    vector<double> m_c;  // ÉÏ¶Ô½ÇÏß
    vector<double> m_p;  // ·Ö½âÖÐ¼ä±äÁ¿
    vector<double> m_q;  // ·Ö½âÖÐ¼ä±äÁ¿
    bool decomposed;
    void decompose();
public:
    // ÌØÊâ¹¹Ôìº¯ÊýÓÃÓÚÈý¶Ô½Ç¾ØÕó
    ChaseMethod(const vector<double>& a, const vector<double>& b,
        const vector<double>& c, const vector<double>& d);
    void solve() override;
};

// Cholesky·Ö½â·¨
class CholeskyDecomposition : public DirectMethod {
private:
    vector<vector<double>> m_L;   // ÏÂÈý½Ç¾ØÕó
    vector<vector<double>> m_LT;  // ×ªÖÃ¾ØÕó
    void choleskyDecompose();
public:
    CholeskyDecomposition(const vector<vector<double>>& A, const vector<double>& b)
        : DirectMethod(A, b) {
        m_L.resize(n, vector<double>(n, 0.0));
        m_LT.resize(n, vector<double>(n, 0.0));
    }
    void solve() override;
};