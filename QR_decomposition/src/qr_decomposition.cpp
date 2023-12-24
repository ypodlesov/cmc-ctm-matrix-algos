#include "qr_decomposition.h"

#include "common.h"

bool QRDecomposition(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const size_t n = A.Size1;
    Q = std::move(A);
    R = TMatrix<double>(A.Size2, A.Size2);
    R.Nullify();
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < j; ++i) {
            double innerProd = TMatrix<double>::InnerProd(Q, i, Q, j);
            R.Data[j * n + i] = innerProd;
            for (size_t k = 0; k < n; ++k) {
                Q.Data[j * n + k] -= Q.Data[i * n + k] * innerProd;
            }
        } 
        double rjj = Q.ColumnNorm2(j);
        R.Data[j * n + j] = rjj;
        for (size_t k = 0; k < n; ++k) {
            Q.Data[j * n + k] /= rjj;
        }
    }
    return true;
}

bool QRDecompositionBlockOptimized(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const size_t blockSize = 32;
    const size_t n = A.Size1;
    Q = std::move(A);
    R = TMatrix<double>(n, n);
    R.Nullify();
    for (size_t j = 0; j < n; ++j) {
        double* jthRColumn = R.GetColumn(j);
        double* jthQColumn = Q.GetColumn(j);
        for (size_t i = 0; i < j; i += blockSize) {
            size_t actualBlockSize = std::min(blockSize, j - i);
            double innerProdVector[actualBlockSize];
            for (size_t k = 0; k < actualBlockSize; ++k) {
                innerProdVector[k] = 0;
            }
            double* ithQColumn = Q.GetColumn(i);
            for (size_t row = 0; row < n; ++row) {
                double coeff = jthQColumn[row];
                for (size_t k = 0; k < actualBlockSize; ++k) {
                    innerProdVector[k] += coeff * ithQColumn[k * n + row];
                }
            }
            for (size_t k = 0; k < actualBlockSize; ++k) {
                jthRColumn[i + k] = innerProdVector[k];
            }

            for (size_t l = 0; l < actualBlockSize; ++l) {
                for (size_t k = 0; k < n; ++k) {
                    jthQColumn[k] -= ithQColumn[l * n + k] * innerProdVector[l];
                }
            }
        }
        double rjj = Q.ColumnNorm2(j);
        R.Data[j * n + j] = rjj;
        for (size_t k = 0; k < n; ++k) {
            Q.Data[j * n + k] /= rjj;
        }
    }
    return true;
}