#include "linear_equations_iterative.h"
iterative_method::iterative_method(const vector<vector<double>>& A, const vector<double> &b, double diff)
    : m_A(A), m_B(b), n(static_cast<int>(b.size())), differ(diff) {}

void iterative_method::print() const {
	for (int i = 0; i < n; ++i) {
		cout << "x" << i << ":  " << ans[i] << endl;
	}
}

Jacobi::Jacobi(const vector<vector<double>>& A, const vector<double>& b, double diff):iterative_method(A,b,diff){}
void Jacobi::solve() {
    ans.resize(n, 0);
    vector<double> prev_ans(n, 0);  // Store the previous iteration result
    const int max_iterations = 1000;
    for (int iter = 0; iter < max_iterations; ++iter) {
        vector<double> new_ans(n, 0);  // Store the current iteration result
        // Jacobi iteration
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += m_A[i][j] * prev_ans[j];
                }
            }
            new_ans[i] = (m_B[i] - sum) / m_A[i][i];
        }
        // Compute the maximum difference between iterations
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(new_ans[i] - prev_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        ans = new_ans;  // Update the final result
        prev_ans = new_ans;  // Prepare for the next iteration
        // Convergence check
        if (max_diff < differ) {
            cout << "Jacobi converged at iteration " << iter + 1 << " iterations" << endl;
            return;
        }
    }
    cout << "Warning: Jacobi did not converge within " << max_iterations << " iterationsiterations" << endl;
}

Gauss_Seidel::Gauss_Seidel(const vector<vector<double>>& A, const vector<double>& b, double diff) :iterative_method(A, b, diff) {};
void Gauss_Seidel::solve() {
    ans.resize(n, 0);
    vector<double> prev_ans(n, 0);
    int max_iterations = 1000;  // Maximum number of iterations

    for (int iter = 0; iter < max_iterations; ++iter) {
        prev_ans = ans;  // Save the previous iteration result
        // Gauss-Seidel iteration
        for (int i = 0; i < n; ++i) {
            double sum = m_B[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= m_A[i][j] * ans[j];
                }
            }
            ans[i] = sum / m_A[i][i];
        }
        // Check convergence
        double max_error = 0.0;
        for (int i = 0; i < n; ++i) {
            max_error = std::max(max_error, std::abs(ans[i] - prev_ans[i]));
        }
        if (max_error < differ) {
            cout << "Gauss-Seidel converged at iteration " << iter + 1 << " iterations" << endl;
            return;
        }
    }
    cout << "Warning: Gauss-Seidel did not converge within the iteration limit" << endl;
}

SOR::SOR(const vector<vector<double>>& A, const vector<double>& b, double diff, double w)
    : iterative_method(A, b, diff), omega(w) {
    // Check whether the relaxation factor is valid
    if (omega <= 0 || omega >= 2) {
        cout << "Warning: omega must be in (0, 2); current value: " << omega << endl;
    }
}

void SOR::solve() {
    vector<double> old_ans(n, 0);  // Store the previous iteration result
    ans.resize(n, 0);

    int iter_count = 0;
    const int max_iterations = 1000;  // Maximum number of iterations

    while (true) {
        old_ans = ans;  // Save the current iteration result

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;

            // Use the latest values
            for (int j = 0; j < i; ++j) {
                sum += m_A[i][j] * ans[j];  // Use values updated in the current iteration
            }
            for (int j = i + 1; j < n; ++j) {
                sum += m_A[i][j] * old_ans[j];  // Use values from the previous iteration
            }

            // SOR iteration formula
            double gauss_seidel = (m_B[i] - sum) / m_A[i][i];
            ans[i] = (1 - omega) * old_ans[i] + omega * gauss_seidel;
        }

        // Compute the maximum difference between iterations
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double diff = abs(ans[i] - old_ans[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        // Check convergence
        if (max_diff < differ) {
            cout << "SOR converged at iteration " << iter_count + 1 << " iterations" << endl;
            break;
        }

        // Check the maximum iteration count
        if (++iter_count >= max_iterations) {
            cout << "Warning: SOR did not converge within " << max_iterations << " iterations" << endl;
            break;
        }
    }
}