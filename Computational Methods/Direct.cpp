#include "Direct.h"

// 高斯消去法实现
void GaussianElimination::solve() {
    vector<vector<double>> matrix(n, vector<double>(n + 1));

    // 构造增广矩阵
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = m_A[i][j];
        }
        matrix[i][n] = m_B[i];
    }

    // 前向消元
    for (int i = 0; i < n; i++) {
        // 列主元选取
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // 行交换
        if (maxRow != i) {
            swap(matrix[i], matrix[maxRow]);
        }

        // 检查主元
        if (abs(matrix[i][i]) < 1e-10) {
            throw runtime_error("矩阵奇异，无法求解");
        }

        // 消元
        for (int k = i + 1; k < n; k++) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j <= n; j++) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }

    // 回代求解
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            ans[i] -= matrix[i][j] * ans[j];
        }
        ans[i] /= matrix[i][i];
    }
}

// LU分解实现
void LU_Decomposition::doolittleDecompose() {
    for (int i = 0; i < n; i++) {
        m_L[i][i] = 1.0;

        // 计算U的第i行
        for (int j = i; j < n; j++) {
            m_U[i][j] = m_A[i][j];
            for (int k = 0; k < i; k++) {
                m_U[i][j] -= m_L[i][k] * m_U[k][j];
            }
        }

        // 计算L的第i列
        for (int j = i + 1; j < n; j++) {
            m_L[j][i] = m_A[j][i];
            for (int k = 0; k < i; k++) {
                m_L[j][i] -= m_L[j][k] * m_U[k][i];
            }
            m_L[j][i] /= m_U[i][i];
        }
    }
}

void LU_Decomposition::solve() {
    doolittleDecompose();

    // 解 L * y = B
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = m_B[i];
        for (int j = 0; j < i; j++) {
            y[i] -= m_L[i][j] * y[j];
        }
        y[i] /= m_L[i][i];
    }

    // 解 U * x = y
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            ans[i] -= m_U[i][j] * ans[j];
        }
        ans[i] /= m_U[i][i];
    }
}

// 追赶法实现
ChaseMethod::ChaseMethod(const vector<double>& a, const vector<double>& b,
    const vector<double>& c, const vector<double>& d)
    : m_a(a), m_c(c), DirectMethod(vector<vector<double>>(), d) {
    // 构造三对角矩阵（仅用于基类初始化）
    m_A.resize(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
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

// Cholesky分解实现
void CholeskyDecomposition::choleskyDecompose() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += m_L[j][k] * m_L[j][k];
                }
                if (m_A[j][j] - sum <= 0) {
                    throw runtime_error("矩阵不正定，无法进行Cholesky分解");
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

    // 计算转置
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m_LT[i][j] = m_L[j][i];
        }
    }
}

void CholeskyDecomposition::solve() {
    choleskyDecompose();

    // 前向替换：L * y = B
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        y[i] = m_B[i];
        for (int j = 0; j < i; j++) {
            y[i] -= m_L[i][j] * y[j];
        }
        y[i] /= m_L[i][i];
    }
    // 后向替换：L^T * x = y
    for (int i = n - 1; i >= 0; i--) {
        ans[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            ans[i] -= m_LT[i][j] * ans[j];
        }
        ans[i] /= m_LT[i][i];
    }
}