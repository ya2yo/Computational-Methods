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
    vector<double> prev_ans(n, 0);  // ´æ´¢ÉÏÒ»´Îµü´ú½á¹û
    const int max_iterations = 1000;
    for (int iter = 0; iter < max_iterations; ++iter) {
        vector<double> new_ans(n, 0);  // ´æ´¢µ±Ç°µü´ú½á¹û
        // Jacobi µü´ú
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += m_A[i][j] * prev_ans[j];
                }
            }
            new_ans[i] = (m_B[i] - sum) / m_A[i][i];
        }
        // ¼ÆËãÁ½´Îµü´ú½á¹ûµÄ×î´ó²îÖµ
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(new_ans[i] - prev_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        ans = new_ans;  // ¸üÐÂ×îÖÕ½á¹û
        prev_ans = new_ans;  // ÎªÏÂÒ»´Îµü´ú×¼±¸
        // ÊÕÁ²ÅÐ¶Ï
        if (max_diff < differ) {
            cout << "Jacobi ·½·¨ÊÕÁ²ÓÚµÚ " << iter + 1 << " ´Îµü´ú" << endl;
            return;
        }
    }
    cout << "¾¯¸æ£ºJacobi ·½·¨Î´ÔÚ " << max_iterations << " ´Îµü´úÄÚÊÕÁ²" << endl;
}

Gauss_Seidel::Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double diff) :iterative_method(A, b, diff) {};
void Gauss_Seidel::solve() {
    ans.resize(n, 0);
    vector<double> prev_ans(n, 0);
    int max_iterations = 1000;  // ×î´óµü´ú´ÎÊý

    for (int iter = 0; iter < max_iterations; ++iter) {
        prev_ans = ans;  // ±£´æÉÏÒ»´Îµü´ú½á¹û
        // Gauss-Seidel µü´ú
        for (int i = 0; i < n; ++i) {
            double sum = m_B[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= m_A[i][j] * ans[j];
                }
            }
            ans[i] = sum / m_A[i][i];
        }
        // ¼ì²éÊÕÁ²
        double max_error = 0.0;
        for (int i = 0; i < n; ++i) {
            max_error = std::max(max_error, std::abs(ans[i] - prev_ans[i]));
        }
        if (max_error < differ) {
            cout << "Gauss-Seidel ÊÕÁ²ÓÚµÚ " << iter + 1 << " ´Îµü´ú" << endl;
            return;
        }
    }
    cout << "¾¯¸æ£ºÎ´ÔÚ×î´óµü´ú´ÎÊýÄÚÊÕÁ²" << endl;
}

SOR::SOR(const vector<vector<double>>& A, const vector<double>& b, double diff, double w)
    : iterative_method(A, b, diff), omega(w) {
    // ¼ì²éËÉ³ÚÒò×ÓÊÇ·ñºÏÀí
    if (omega <= 0 || omega >= 2) {
        cout << "¾¯¸æ£ºËÉ³ÚÒò×Ó omega Ó¦ÔÚ (0,2) ·¶Î§ÄÚ£¬µ±Ç°Öµ: " << omega << endl;
    }
}

void SOR::solve() {
    vector<double> old_ans(n, 0);  // ´æ´¢ÉÏÒ»´Îµü´ú½á¹û
    ans.resize(n, 0);

    int iter_count = 0;
    const int max_iterations = 1000;  // ×î´óµü´ú´ÎÊý

    while (true) {
        old_ans = ans;  // ±£´æµ±Ç°½á¹û

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;

            // Ê¹ÓÃ×îÐÂÖµ¼ÆËã
            for (int j = 0; j < i; ++j) {
                sum += m_A[i][j] * ans[j];  // Ê¹ÓÃµ±Ç°µü´úÒÑ¸üÐÂµÄÖµ
            }
            for (int j = i + 1; j < n; ++j) {
                sum += m_A[i][j] * old_ans[j];  // Ê¹ÓÃÉÏÒ»´Îµü´úµÄÖµ
            }

            // SOR µü´ú¹«Ê½: x_i^(k+1) = (1-¦Ø)x_i^(k) + ¦Ø(b_i - ¦²)/a_ii
            double gauss_seidel = (m_B[i] - sum) / m_A[i][i];
            ans[i] = (1 - omega) * old_ans[i] + omega * gauss_seidel;
        }

        // ¼ÆËãÁ½´Îµü´ú½á¹ûµÄ×î´ó²îÖµ
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(ans[i] - old_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        // ¼ì²éÊÕÁ²
        if (max_diff < differ) {
            cout << "SOR ·½·¨ÊÕÁ²ÓÚµÚ " << iter_count + 1 << " ´Îµü´ú" << endl;
            break;
        }

        // ¼ì²é×î´óµü´ú´ÎÊý
        if (++iter_count >= max_iterations) {
            cout << "¾¯¸æ£ºSOR ·½·¨ÔÚ " << max_iterations << " ´Îµü´úÄÚÎ´ÊÕÁ²" << endl;
            break;
        }
    }
}