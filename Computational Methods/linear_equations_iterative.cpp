#include "linear_equations_iterative.h"
iterative_method::iterative_method(const vector<vector<double>>& A, const vector<double> &b,double diff):m_A(A),m_B(b),n(b.size()),differ(diff) {}

void iterative_method::print() const {
	for (size_t i = 0; i < n; ++i) {
		cout << "x" << i << ":  " << ans[i] << endl;
	}
}

Jacobi::Jacobi(const vector<vector<double>>& A, const vector<double>& b, double diff):iterative_method(A,b,diff){}
void Jacobi::solve() {
    ans.resize(n, 0);
    vector<double> prev_ans(n, 0);  // 存储上一次迭代结果
    const int max_iterations = 1000;
    for (int iter = 0; iter < max_iterations; ++iter) {
        vector<double> new_ans(n, 0);  // 存储当前迭代结果
        // Jacobi 迭代
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += m_A[i][j] * prev_ans[j];
                }
            }
            new_ans[i] = (m_B[i] - sum) / m_A[i][i];
        }
        // 计算两次迭代结果的最大差值
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(new_ans[i] - prev_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        ans = new_ans;  // 更新最终结果
        prev_ans = new_ans;  // 为下一次迭代准备
        // 收敛判断
        if (max_diff < differ) {
            cout << "Jacobi 方法收敛于第 " << iter + 1 << " 次迭代" << endl;
            return;
        }
    }
    cout << "警告：Jacobi 方法未在 " << max_iterations << " 次迭代内收敛" << endl;
}

Gauss_Seidel::Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double diff) :iterative_method(A, b, diff) {};
void Gauss_Seidel::solve() {
    ans.resize(n, 0);
    vector<double> prev_ans(n, 0);
    int max_iterations = 1000;  // 最大迭代次数

    for (int iter = 0; iter < max_iterations; ++iter) {
        prev_ans = ans;  // 保存上一次迭代结果
        // Gauss-Seidel 迭代
        for (int i = 0; i < n; ++i) {
            double sum = m_B[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= m_A[i][j] * ans[j];
                }
            }
            ans[i] = sum / m_A[i][i];
        }
        // 检查收敛
        double max_error = 0.0;
        for (int i = 0; i < n; ++i) {
            max_error = std::max(max_error, std::abs(ans[i] - prev_ans[i]));
        }
        if (max_error < differ) {
            cout << "Gauss-Seidel 收敛于第 " << iter + 1 << " 次迭代" << endl;
            return;
        }
    }
    cout << "警告：未在最大迭代次数内收敛" << endl;
}

SOR::SOR(const vector<vector<double>>& A, const vector<double>& b, double diff, double w)
    : iterative_method(A, b, diff), omega(w) {
    // 检查松弛因子是否合理
    if (omega <= 0 || omega >= 2) {
        cout << "警告：松弛因子 omega 应在 (0,2) 范围内，当前值: " << omega << endl;
    }
}

void SOR::solve() {
    vector<double> old_ans(n, 0);  // 存储上一次迭代结果
    ans.resize(n, 0);

    int iter_count = 0;
    const int max_iterations = 1000;  // 最大迭代次数

    while (true) {
        old_ans = ans;  // 保存当前结果

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;

            // 使用最新值计算
            for (int j = 0; j < i; ++j) {
                sum += m_A[i][j] * ans[j];  // 使用当前迭代已更新的值
            }
            for (int j = i + 1; j < n; ++j) {
                sum += m_A[i][j] * old_ans[j];  // 使用上一次迭代的值
            }

            // SOR 迭代公式: x_i^(k+1) = (1-ω)x_i^(k) + ω(b_i - Σ)/a_ii
            double gauss_seidel = (m_B[i] - sum) / m_A[i][i];
            ans[i] = (1 - omega) * old_ans[i] + omega * gauss_seidel;
        }

        // 计算两次迭代结果的最大差值
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(ans[i] - old_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        // 检查收敛
        if (max_diff < differ) {
            cout << "SOR 方法收敛于第 " << iter_count + 1 << " 次迭代" << endl;
            break;
        }

        // 检查最大迭代次数
        if (++iter_count >= max_iterations) {
            cout << "警告：SOR 方法在 " << max_iterations << " 次迭代内未收敛" << endl;
            break;
        }
    }
}