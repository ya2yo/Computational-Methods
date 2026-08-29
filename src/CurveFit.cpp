#include "CurveFit.h"
double Data::evaluate(double z) const {
    if (factor.empty()) {
        throw runtime_error("ÉÐÎ´¼ÆËãÄâºÏÏµÊý (factor)!");
    }
    double result = 0.0;

    for (size_t i = 0; i < factor.size(); ++i) {
        result += factor[i] * power(z, i);
    }
    return result;
}

void Ordinary_Method::calculateFactor(int n) {
    factor.clear();

    if (n < 0 || static_cast<std::size_t>(n) >= x.size()) {
        throw runtime_error("¶àÏîÊ½´ÎÊý n ²»ºÏÀí (Ó¦ÔÚ 0 µ½ Êý¾Ýµã¸öÊý-1 Ö®¼ä)");
    }

    Matrix A(n + 1, Vector(n + 1, 0.0));
    Vector b(n + 1, 0.0);
    int m = x.size(); 


    for (int i = 0; i <= n; ++i) {
        for (int j = i; j <= n; ++j) { 
            double sum_pow = 0.0;
            for (int k = 0; k < m; ++k) {
                sum_pow += power(x[k], i + j);
            }
            A[i][j] = sum_pow;
            if (i != j) {
                A[j][i] = sum_pow; 
            }
        }

        for (int k = 0; k < m; ++k) {
            b[i] += y[k] * power(x[k], i);
        }
    }


    factor = solveLinearSystem(A, b);
}

void Ordinary_Method::calculate(int n)
{
    calculateFactor(n);
}

void Orthogonal_Method::calculateCoefficients(int n) {
    factor.clear();
    alpha.clear();
    beta.clear();

    if (n < 0 || static_cast<std::size_t>(n) >= x.size()) {
        throw runtime_error("¶àÏîÊ½´ÎÊý n ²»ºÏÀí (Ó¦ÔÚ 0 µ½ Êý¾Ýµã¸öÊý-1 Ö®¼ä)");
    }

    int m = x.size(); 

    factor.resize(n + 1);
    alpha.resize(n + 1, 0.0); 
    beta.resize(n + 1, 0.0); 


    Vector phi_curr(m);
    Vector phi_prev(m, 0.0); 
    Vector phi_pprev(m, 0.0);

    double sum_y_phi = 0.0;
    double sum_phi_curr_sq = 0.0;
    double sum_x_phi_sq = 0.0;

    double sum_phi_prev_sq = 0.0; 
    double sum_phi_pprev_sq = 0.0; 

    sum_y_phi = 0.0;
    sum_x_phi_sq = 0.0;
    sum_phi_curr_sq = (double)m;

    for (int k = 0; k < m; ++k) {
        phi_curr[k] = 1.0;
        sum_y_phi += y[k];     
        sum_x_phi_sq += x[k];  
    }


    factor[0] = sum_y_phi / sum_phi_curr_sq;

    alpha[0] = sum_x_phi_sq / sum_phi_curr_sq;


    phi_prev = phi_curr;
    sum_phi_prev_sq = sum_phi_curr_sq;

    for (int i = 1; i <= n; ++i) {


        if (i >= 2) {
            if (sum_phi_pprev_sq < 1e-12) {
                throw runtime_error("Õý½»¶àÏîÊ½ÄâºÏÖÐ Psi_{i-2} Æ½·½ºÍ½Ó½üÁã£¬ÎÞ·¨¼ÆËã beta¡£");
            }
            beta[i - 1] = sum_phi_prev_sq / sum_phi_pprev_sq;
        }

        sum_phi_curr_sq = 0.0;
        sum_y_phi = 0.0;
        sum_x_phi_sq = 0.0;

        for (int k = 0; k < m; ++k) {
            double current_beta = (i == 1) ? 0.0 : beta[i - 1];

            // µÝÍÆ¼ÆËã
            phi_curr[k] = (x[k] - alpha[i - 1]) * phi_prev[k] - current_beta * phi_pprev[k];

            // ¼ÆËãÐÂµÄºÍ
            sum_phi_curr_sq += phi_curr[k] * phi_curr[k];
            sum_y_phi += y[k] * phi_curr[k];
            sum_x_phi_sq += x[k] * phi_curr[k] * phi_curr[k];
        }

        // ¼ì²é·ÖÄ¸
        if (sum_phi_curr_sq < 1e-12) {
            throw runtime_error("Õý½»¶àÏîÊ½ÄâºÏÖÐ Sum [Psi_i]^2 ½Ó½üÁã£¬ÄâºÏÖÕÖ¹¡£");
        }

        // 3. ¼ÆËã c_i
        factor[i] = sum_y_phi / sum_phi_curr_sq;

        // 4. ¼ÆËã alpha_i (Èç¹ûÐèÒªÏÂÒ»ÂÖµÝÍÆ)
        if (i < n) {
            alpha[i] = sum_x_phi_sq / sum_phi_curr_sq;
        }


        phi_pprev = phi_prev;
        phi_prev = phi_curr; 

        sum_phi_pprev_sq = sum_phi_prev_sq; 
        sum_phi_prev_sq = sum_phi_curr_sq;
    }
}

void Orthogonal_Method::calculate(int n) {
    calculateCoefficients(n);
}


double Orthogonal_Method::evaluate(double z) const {
    if (factor.empty()) {
        throw runtime_error("ÉÐÎ´¼ÆËãÄâºÏÏµÊý (factor)!");
    }

    int n = factor.size() - 1;

    double Psi_pprev = 0.0;
    double Psi_prev = 1.0;  
    double Psi_curr = 0.0;

    double result = factor[0] * Psi_prev; 

    for (int i = 1; i <= n; ++i) {
        double current_beta = (i == 1) ? 0.0 : beta[i - 1];

        Psi_curr = (z - alpha[i - 1]) * Psi_prev - current_beta * Psi_pprev;

        result += factor[i] * Psi_curr;

        Psi_pprev = Psi_prev;
        Psi_prev = Psi_curr;
    }

    return result;
}