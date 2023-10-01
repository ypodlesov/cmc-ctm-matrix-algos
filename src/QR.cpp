#include "QR.h"

bool QR_decomposition(matrix<double>& A, matrix<double>& Q, triangular_matrix<double, upper>& R) {
    if (A.size1() < A.size2()) {
        return false;
    }
    const size_t m = A.size1();
    const size_t n = A.size2();
    Q = matrix<double>(A.size1(), A.size2());
    R = triangular_matrix<double, upper>(A.size2(), A.size2());
    for (size_t k = 0; k < n; ++k) {
        const double norm = norm_2(vector<double>(matrix_vector_slice<matrix<double>>(A, slice(0, 1, m), slice(k, 0, m))));
        R(k, k) = norm;
        for (size_t i = 0; rough_lt(0.0, norm) && i < m; ++i) {
            Q(i, k) = A(i, k) / norm;
        }
        for (size_t j = k + 1; j < n; ++j) {
            for (size_t i = 0; i < m; ++i) {
                R(k, j) += Q(i, k) * A(i, j); 
            }
            for (size_t i = 0; i < m; ++i) {
                A(i, j) -= Q(i, k) * R(k, j);
            }
        }
    }
    return true;
} 