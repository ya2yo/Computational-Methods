#include "Direct.h"

#include <cmath>
#include <utility>

using namespace std;

// Gaussian elimination implementation
void GaussianElimination::solve() {
    vector<vector<double>> matrix(n, vector<double>(n + 1));

    // Build the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) matrix[i][j] = m_A[i][j];
        matrix[i][n] = m_B[i];
    }

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(matrix[k][i]) > abs(matrix[maxRow][i])) maxRow = k;
        }
        if (maxRow != i) swap(matrix[i], matrix[maxRow]);
        if (abs(matrix[i][i]) < 1e-10) {
            throw runtime_error("The matrix is singular; no solution can be computed");
        }
        for (int k = i + 1; k < n; ++k) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j <= n; ++j) matrix[k][j] -= factor * matrix[i][j];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        ans[i] = matrix[i][n];
        for (int j = i + 1; j < n; ++j) ans[i] -= matrix[i][j] * ans[j];
        ans[i] /= matrix[i][i];
    }
}

// LU decomposition implementation
void LU_Decomposition::doolittleDecompose() {
    for (int i = 0; i < n; ++i) {
        m_L[i][i] = 1.0;
        for (int j = i; j < n; ++j) {
            m_U[i][j] = m_A[i][j];
            for (int k = 0; k < i; ++k) m_U[i][j] -= m_L[i][k] * m_U[k][j];
        }
        for (int j = i + 1; j < n; ++j) {
            m_L[j][i] = m_A[j][i];
            for (int k = 0; k < i; ++k) m_L[j][i] -= m_L[j][k] * m_U[k][i];
            m_L[j][i] /= m_U[i][i];
        }
    }
}

void LU_Decomposition::solve() {
    doolittleDecompose();
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = m_B[i];
        for (int j = 0; j < i; ++j) y[i] -= m_L[i][j] * y[j];
        y[i] /= m_L[i][i];
    }
    for (int i = n - 1; i >= 0; --i) {
        ans[i] = y[i];
        for (int j = i + 1; j < n; ++j) ans[i] -= m_U[i][j] * ans[j];
        ans[i] /= m_U[i][i];
    }
}

// Thomas algorithm implementation
ChaseMethod::ChaseMethod(const vector<double>& a, const vector<double>& b,
    const vector<double>& c, const vector<double>& d)
    : DirectMethod(vector<vector<double>>(), d), m_a(a), m_c(c) {
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
    for (int i = 1; i < n; ++i) {
        double denom = m_A[i][i] - m_a[i] * m_p[i - 1];
        if (i < n - 1) m_p[i] = m_c[i] / denom;
        m_q[i] = (m_B[i] - m_a[i] * m_q[i - 1]) / denom;
    }
    decomposed = true;
}

void ChaseMethod::solve() {
    if (!decomposed) decompose();
    ans[n - 1] = m_q[n - 1];
    for (int i = n - 2; i >= 0; --i) ans[i] = m_q[i] - m_p[i] * ans[i + 1];
}

// Cholesky decomposition implementation
void CholeskyDecomposition::choleskyDecompose() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            if (j == i) {
                for (int k = 0; k < j; ++k) sum += m_L[j][k] * m_L[j][k];
                if (m_A[j][j] - sum <= 0) {
                    throw runtime_error("The matrix is not positive definite; Cholesky decomposition failed");
                }
                m_L[j][j] = sqrt(m_A[j][j] - sum);
            } else {
                for (int k = 0; k < j; ++k) sum += m_L[i][k] * m_L[j][k];
                m_L[i][j] = (m_A[i][j] - sum) / m_L[j][j];
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) m_LT[i][j] = m_L[j][i];
    }
}

void CholeskyDecomposition::solve() {
    choleskyDecompose();
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = m_B[i];
        for (int j = 0; j < i; ++j) y[i] -= m_L[i][j] * y[j];
        y[i] /= m_L[i][i];
    }
    for (int i = n - 1; i >= 0; --i) {
        ans[i] = y[i];
        for (int j = i + 1; j < n; ++j) ans[i] -= m_LT[i][j] * ans[j];
        ans[i] /= m_LT[i][i];
    }
}
