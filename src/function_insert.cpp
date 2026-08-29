#include "function_insert.h"



// ¼ÆËãÀ­¸ñÀÊÈÕ»ùº¯Êý L_k(x)
double lagrangeBasis(int k, double val, const vector<double>& x_nodes) {
    double result = 1.0;
    for (std::size_t i = 0; i < x_nodes.size(); ++i) {
        if (i != static_cast<std::size_t>(k)) {
            // L_k(x) = Prod_{i!=k} (x - x_i) / (x_k - x_i)
            result *= (val - x_nodes[i]) / (x_nodes[k] - x_nodes[i]);
        }
    }
    return result;
}

double LagrangeInterpolation::interpolate(double val) {
    double P_n = 0.0;
    for (int k = 0; k < n; k++) {
        // P_n(x) = Sum_{k=0}^n y_k * L_k(x)
        P_n += y[k] * lagrangeBasis(k, val, x);
    }
    return P_n;
}


void NewtonInterpolation::computeDividedDiffs() {
    
    divided_diffs = y; 

    for (int k = 1; k < n; k++) { 
        for (int i = n - 1; i >= k; i--) { 
            
            divided_diffs[i] = (divided_diffs[i] - divided_diffs[i - 1]) / (x[i] - x[i - k]);
        }
    }
}

double NewtonInterpolation::interpolate(double val) {
    double P_n = divided_diffs[0];
    double product = 1.0;

    for (int k = 1; k < n; k++) {
        
        product *= (val - x[k - 1]);
        P_n += divided_diffs[k] * product;
    }
    return P_n;
}



void HermiteInterpolation::computeHermiteCoeffs() {
    
    vector<double> z(2 * n);
    vector<double> f(2 * n); 

    for (int i = 0; i < n; i++) {
        z[2 * i] = x[i];
        z[2 * i + 1] = x[i];
        f[2 * i] = y[i];
        f[2 * i + 1] = y[i];
    }

    
    H_coeffs.resize(2 * n);
    for (int i = 0; i < 2 * n; i++) {
        H_coeffs[i] = f[i];
    }

    
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            H_coeffs[2 * i + 1] = y_prime[i];
        }
        else {
            
            H_coeffs[2 * i] = (H_coeffs[2 * i] - H_coeffs[2 * i - 1]) / (z[2 * i] - z[2 * i - 1]);
            H_coeffs[2 * i + 1] = y_prime[i]; // f[x_i, x_i] = y'_i
        }
    }

    
    for (int k = 2; k < 2 * n; k++) { 
        for (int i = 2 * n - 1; i >= k; i--) {
            
            H_coeffs[i] = (H_coeffs[i] - H_coeffs[i - 1]) / (z[i] - z[i - k]);
        }
    }
    
}

double HermiteInterpolation::interpolate(double val) {
    
    vector<double> z(2 * n);
    for (int i = 0; i < n; i++) {
        z[2 * i] = x[i];
        z[2 * i + 1] = x[i];
    }

    double H_val = H_coeffs[0]; 
    double product = 1.0;

    for (int k = 1; k < 2 * n; k++) {
        
        product *= (val - z[k - 1]);
        H_val += H_coeffs[k] * product;
    }
    return H_val;
}



void CubicSplineInterpolation::solveTridiagonalSystem() {
    if (n < 3) return; 


    
    int dim = n - 2; 
    if (dim <= 0) return;

    
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    
    vector<double> a(dim);
    
    vector<double> b(dim);
    
    vector<double> c(dim);
    
    vector<double> d(dim);

    for (int i = 0; i < dim; i++) {
        int k = i + 1; 

        
        if (k > 0) a[i] = h[k - 1]; 
        b[i] = 2.0 * (h[k - 1] + h[k]); 
        if (k < n - 2) c[i] = h[k]; 

        
        double delta_y1 = (y[k + 1] - y[k]) / h[k];
        double delta_y2 = (y[k] - y[k - 1]) / h[k - 1];
        d[i] = 6.0 * (delta_y1 - delta_y2);
    }
    vector<double> P(dim);
    vector<double> Q(dim);

    
    P[0] = c[0] / b[0];
    Q[0] = d[0] / b[0];

    for (int i = 1; i < dim; i++) {
        double denom = b[i] - a[i] * P[i - 1];
        if (abs(denom) < 1e-10) {
            throw runtime_error("Èý´ÎÑùÌõ²åÖµ¾ØÕóÆæÒì£¬ÎÞ·¨Çó½â");
        }
        if (i < dim - 1) {
            P[i] = c[i] / denom;
        }
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom;
    }

    
    vector<double> temp_M(dim);
    temp_M[dim - 1] = Q[dim - 1]; // M_{n-2}

    for (int i = dim - 2; i >= 0; i--) {
        temp_M[i] = Q[i] - P[i] * temp_M[i + 1]; // M_{i+1}
    }

    for (int i = 0; i < dim; i++) {
        M[i + 1] = temp_M[i];
    }
}


int CubicSplineInterpolation::findInterval(double val) {
    if (val < x[0] || val > x[n - 1]) {
        throw out_of_range("²åÖµµã³¬³ö¸ø¶¨½Úµã·¶Î§");
    }

    for (int i = 0; i < n - 1; i++) {
        if (val >= x[i] && val <= x[i + 1]) {
            return i;
        }
    }
    return n - 2;
}

double CubicSplineInterpolation::interpolate(double val) {
    if (n < 3) {
        if (n == 2) {
            return y[0] + (y[1] - y[0]) * (val - x[0]) / (x[1] - x[0]);
        }
        else if (n == 1) {
            return y[0];
        }
        throw invalid_argument("½ÚµãÊýÁ¿²»×ã");
    }
    int i = findInterval(val);
    int i_plus_1 = i + 1;

    double h_i = x[i_plus_1] - x[i];
    double u = (x[i_plus_1] - val) / h_i;
    double v = (val - x[i]) / h_i;

    
    double S_val =
        M[i] * u * (u * u - 1.0) * (h_i * h_i / 6.0) +
        M[i_plus_1] * v * (v * v - 1.0) * (h_i * h_i / 6.0) +
        y[i] * u +
        y[i_plus_1] * v;

    return S_val;
}