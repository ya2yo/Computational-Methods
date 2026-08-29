#include "Direct.h"

#include <cmath>
#include <utility>

// ¸ßË¹ÏûÈ¥·¨ÊµÏÖ
void GaussianElimination::solve() {
    vector<vector<double>> matrix(n, vector<double>(n + 1));

    // ¹¹ÔìÔö¹ã¾ØÕó
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = m_A[i][j];
        }
        matrix[i][n] = m_B[i];
    }

    // Ç°ÏòÏûÔª
    for (int i = 0; i < n; ++i) {
        // ÁÐÖ÷ÔªÑ¡È¡
        int maxRow = static_cast<int>(i);
        for (int k = i + 1; k < n; ++k) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // ÐÐ½»»»
        if (maxRow != i) {
            swap(matrix[i], matrix[maxRow]);
        }

        // ¼ì²éÖ÷Ôª
        if (abs(matrix[i][i]) < 1e-10) {
            throw runtime_error("¾ØÕóÆæÒì£¬ÎÞ·¨Çó½â");
        }

        // ÏûÔª
        for (int k = i + 1; k < n; ++k) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = static_cast<int>(i); j <= n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }

    // »Ø´úÇó½â
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = matrix[i][n];
        for (int j = static_cast<int>(i) + 1; j < n; ++j) {
            ans[i] -= matrix[i][j] * ans[j];
        }
        ans[i] /= matrix[i][i];
    }
}

// LU·Ö½âÊµÏÖ
void LU_Decomposition::doolittleDecompose() {
    for (int i = 0; i < n; ++i) {
        m_L[i][i] = 1.0;

        // ¼ÆËãUµÄµÚiÐÐ
        for (int j = static_cast<int>(i); j < n; ++j) {
            m_U[i][j] = m_A[i][j];
            for (int k = 0; k < i; ++k) {
                m_U[i][j] -= m_L[i][k] * m_U[k][j];
            }
        }

        // ¼ÆËãLµÄµÚiÁÐ
        for (int j = static_cast<int>(i) + 1; j < n; ++j) {
            m_L[j][i] = m_A[j][i];
            for (int k = 0; k < i; ++k) {
                m_L[j][i] -= m_L[j][k] * m_U[k][i];
            }
            m_L[j][i] /= m_U[i][i];
        }
    }
}

void LU_Decomposition::solve() {
    doolittleDecompose();

    // ½â L * y = B
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = m_B[i];
        for (int j = 0; j < static_cast<int>(i); ++j) {
            y[i] -= m_L[i][j] * y[j];
        }
        y[i] /= m_L[i][i];
    }

    // ½â U * x = y
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = y[i];
        for (int j = static_cast<int>(i) + 1; j < n; ++j) {
            ans[i] -= m_U[i][j] * ans[j];
        }
        ans[i] /= m_U[i][i];
    }
}

// ×·¸Ï·¨ÊµÏÖ
ChaseMethod::ChaseMethod(const vector<double>& a, const vector<double>& b,
    const vector<double>& c, const vector<double>& d)
    : DirectMethod(vector<vector<double>>(), d), m_a(a), m_c(c) {
    // ¹¹ÔìÈý¶Ô½Ç¾ØÕó£¨½öÓÃÓÚ»ùÀà³õÊ¼»¯£©
    m_A.resize(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        if (i > 0) m_A[i][i - 1] = a[i];
        m_A[i][i] = b[i];
        if (i < n - 1) m_A[i][i + 1] = c[i];
    }

    m_p.resize(n);
    m_q.resize(n);
    decomposed = false;
}

void ChaseMethod::decompose() {
    m_p[0] = m_c[0] / m_A[0][0];
    m_q[0] = m_B[0] / m_A[0][0];

    for (int i = 1; i < n; i++) {
        double denom = m_A[i][i] - m_a[i] * m_p[i - 1];
        if (i < n - 1) {
            m_p[i] = m_c[i] / denom;
        }
        m_q[i] = (m_B[i] - m_a[i] * m_q[i - 1]) / denom;
    }

    decomposed = true;
}

void ChaseMethod::solve() {
    if (!decomposed) {
        decompose();
    }

    ans[n - 1] = m_q[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        ans[i] = m_q[i] - m_p[i] * ans[i + 1];
    }
}

// Cholesky·Ö½âÊµÏÖ
void CholeskyDecomposition::choleskyDecompose() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += m_L[j][k] * m_L[j][k];
                }
                if (m_A[j][j] - sum <= 0) {
                    throw runtime_error("¾ØÕó²»Õý¶¨£¬ÎÞ·¨½øÐÐCholesky·Ö½â");
                }
                m_L[j][j] = sqrt(m_A[j][j] - sum);
            }
            else {
                for (int k = 0; k < j; k++) {
                    sum += m_L[i][k] * m_L[j][k];
                }
                m_L[i][j] = (m_A[i][j] - sum) / m_L[j][j];
            }
        }
    }

    // ¼ÆËã×ªÖÃ
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            m_LT[i][j] = m_L[j][i];
        }
    }
}

void CholeskyDecomposition::solve() {
    choleskyDecompose();

    // Ç°ÏòÌæ»»£ºL * y = B
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = m_B[i];
        for (int j = 0; j < static_cast<int>(i); ++j) {
            y[i] -= m_L[i][j] * y[j];
        }
        y[i] /= m_L[i][i];
    }
    // ºóÏòÌæ»»£ºL^T * x = y
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = y[i];
        for (int j = static_cast<int>(i) + 1; j < n; ++j) {
            ans[i] -= m_LT[i][j] * ans[j];
        }
        ans[i] /= m_LT[i][i];
    }
}