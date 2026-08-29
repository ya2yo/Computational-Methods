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
            throw invalid_argument("²åÖµ½ÚµãÊýÁ¿²»Æ¥Åä»òÎª¿Õ");
        }
    }

    virtual ~FunctionInsert() {}


    virtual double interpolate(double val) = 0;


    void printNodes() const {
        cout << "²åÖµ½Úµã£º" << endl;
        for (int i = 0; i < n; i++) {
            cout << "(" << x[i] << ", " << y[i] << ") ";
        }
        cout << endl;
    }
};

// --- À­¸ñÀÊÈÕ²åÖµ ---
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
            throw invalid_argument("ºÕÃ×ÌØ²åÖµÒªÇó½Úµã¡¢º¯ÊýÖµºÍµ¼ÊýÖµÊýÁ¿Ò»ÖÂ");
        }
        computeHermiteCoeffs();
    }


    double interpolate(double val) override;
};

//Èý´ÎÑùÌõº¯Êý²åÖµ (²ÉÓÃ±ß½çÌõ¼þ M0=Mn=0)
class CubicSplineInterpolation : public FunctionInsert {
private:
    vector<double> M; 

    void solveTridiagonalSystem();

public:
    CubicSplineInterpolation(const vector<double>& x_nodes, const vector<double>& y_nodes)
        : FunctionInsert(x_nodes, y_nodes) {
        M.resize(n, 0.0);
        if (n < 3) {
            throw invalid_argument("Èý´ÎÑùÌõ²åÖµÖÁÉÙÐèÒª3¸ö½Úµã");
        }
        solveTridiagonalSystem(); 
    }
    double interpolate(double val) override;
    int findInterval(double val);
};