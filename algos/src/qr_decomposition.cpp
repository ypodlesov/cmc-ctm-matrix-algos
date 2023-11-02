#include "qr_decomposition.h"

#include "common.h"

bool QRDecomposition(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.GetSize1() < A.GetSize2()) {
        return false;
    }
    const size_t m = A.GetSize1();
    const size_t n = A.GetSize2();
    Q = TMatrix<double>(A.GetSize1(), A.GetSize2());
    Q.Nullify();
    R = TMatrix<double>(A.GetSize2(), A.GetSize2());
    R.Nullify();
    for (size_t k = 0; k < n; ++k) {
        const double norm = TMatrix<double>::ColumnNorm2(A, k);
        R(k, k) = norm;
        for (size_t i = 0; !RoughEq(0.0, norm) && i < m; ++i) {
            Q(i, k) = A(i, k) / norm;
        }
        for (size_t j = k + 1; j < n; ++j) {
            R(k, j) += TMatrix<double>::InnerProd(Q, k, A, j);
            for (size_t i = 0; i < m; ++i) {
                A(i, j) -= Q(i, k) * R(k, j);
            }
        }
    }
    return true;
}

// bool QRDecompositionBlockOptimized(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
//     if (A.GetSize1() < A.GetSize2()) {
//         return false;
//     }
//     const size_t m = A.GetSize1();
//     const size_t n = A.GetSize2();
//     Q = TMatrix<double>(A.GetSize1(), A.GetSize2());
//     Q.Nullify();
//     R = TMatrix<double>(A.GetSize2(), A.GetSize2());
//     R.Nullify();
//     for (size_t k = 0; k < n; ++k) {
//         const double norm = TMatrix<double>::ColumnNorm2(A, k);
//         R(k, k) = norm;
//         for (size_t i = 0; !RoughEq(0.0, norm) && i < m; ++i) {
//             Q(i, k) = A(i, k) / norm;
//         }
//         for (size_t j = k + 1; j < n; ++j) {
//             auto v = TVector<double>(Q, k);
//             auto u = TVector<double>(A, j);
//             double rkj = 0.0;
//             rkj += TVector<double>::InnerProd(v, u);
//             for (size_t i = 0; i < m; ++i) {
//                 u(i) -= v(i) * rkj;
//             }
//             R(k, j) = rkj;
//             A.AssignColumn(j, u);
//         }
//     }
//     return true;
// }