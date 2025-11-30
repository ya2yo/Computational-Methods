#include "CurveFit.h"
double Data::evaluate(double z) const {
    if (factor.empty()) {
        throw runtime_error("尚未计算拟合系数 (factor)!");
    }
    double result = 0.0;

    for (size_t i = 0; i < factor.size(); ++i) {
        result += factor[i] * power(z, i);
    }
    return result;
}

void Ordinary_Method::calculateFactor(int n) {
    factor.clear();

    if (n < 0 || n >= x.size()) {
        throw runtime_error("多项式次数 n 不合理 (应在 0 到 数据点个数-1 之间)");
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

    if (n < 0 || n >= x.size()) {
        throw runtime_error("多项式次数 n 不合理 (应在 0 到 数据点个数-1 之间)");
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
                throw runtime_error("正交多项式拟合中 Psi_{i-2} 平方和接近零，无法计算 beta。");
            }
            beta[i - 1] = sum_phi_prev_sq / sum_phi_pprev_sq;
        }

        sum_phi_curr_sq = 0.0;
        sum_y_phi = 0.0;
        sum_x_phi_sq = 0.0;

        for (int k = 0; k < m; ++k) {
            double current_beta = (i == 1) ? 0.0 : beta[i - 1];

            // 递推计算
            phi_curr[k] = (x[k] - alpha[i - 1]) * phi_prev[k] - current_beta * phi_pprev[k];

            // 计算新的和
            sum_phi_curr_sq += phi_curr[k] * phi_curr[k];
            sum_y_phi += y[k] * phi_curr[k];
            sum_x_phi_sq += x[k] * phi_curr[k] * phi_curr[k];
        }

        // 检查分母
        if (sum_phi_curr_sq < 1e-12) {
            throw runtime_error("正交多项式拟合中 Sum [Psi_i]^2 接近零，拟合终止。");
        }

        // 3. 计算 c_i
        factor[i] = sum_y_phi / sum_phi_curr_sq;

        // 4. 计算 alpha_i (如果需要下一轮递推)
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
        throw runtime_error("尚未计算拟合系数 (factor)!");
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