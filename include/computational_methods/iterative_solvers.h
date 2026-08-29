#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

// Iterative methods for linear systems.
class iterative_method {
protected:
    vector<vector<double>> m_A;
    vector<double> m_B;
    vector<double> ans;
    int n;
    double differ;
public:
    iterative_method(const vector<vector<double>>& A, const vector<double>& b, double tolerance);
    virtual void solve() = 0;
    void print() const;
    vector<double> getSolution() const { return ans; }
    virtual ~iterative_method() {}
};

class Jacobi : public iterative_method {
public:
    Jacobi(const vector<vector<double>>& A, const vector<double>& b, double tolerance);
    void solve() override;
};

class Gauss_Seidel : public iterative_method {
public:
    Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double tolerance);
    void solve() override;
};

// Successive over-relaxation iteration.
class SOR : public iterative_method {
private:
    double omega;
public:
    SOR(const vector<vector<double>>& A, const vector<double>& b, double tolerance, double relaxation);
    void solve() override;
};
