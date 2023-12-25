#include "qr_decomposition.h"

#include "common.h"

bool QRDecomposition(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const int n = A.Size1;
    Q = std::move(A);
    R = TMatrix<double>(A.Size2, A.Size2);
    R.Nullify();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) {
            double innerProd = TMatrix<double>::InnerProd(Q, i, Q, j);
            R.Data[j * n + i] = innerProd;
            for (int k = 0; k < n; ++k) {
                Q.Data[j * n + k] -= Q.Data[i * n + k] * innerProd;
            }
        } 
        double rjj = Q.ColumnNorm2(j);
        R.Data[j * n + j] = rjj;
        for (int k = 0; k < n; ++k) {
            Q.Data[j * n + k] /= rjj;
        }
    }
    return true;
}

bool QRDecompositionBlas(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const int increment = 1; 
    const int n = A.Size1;
    Q = std::move(A);
    R = TMatrix<double>(A.Size2, A.Size2);
    R.Nullify();
    for (int j = 0; j < n; ++j) {
        double* jthQColumn = Q.GetColumn(j);
        double* jthRColumn = R.GetColumn(j);
        for (int i = 0; i < j; ++i) {
            double* ithQColumn = Q.GetColumn(i);
            double innerProd = ddot_(&n, jthQColumn, &increment, ithQColumn, &increment);
            jthRColumn[i] = innerProd;
            innerProd *= -1;
            daxpy_(&n, &innerProd, ithQColumn, &increment, jthQColumn, &increment);
        } 
        double rjj = dnrm2_(&n, jthQColumn, &increment);
        jthRColumn[j] = rjj;
        const double coeff = 1 / rjj;
        dscal_(&n, &coeff, jthQColumn, &increment);
    }
    return true;
}

bool QRDecompositionBlockOptimized(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const int blockSize = 8;
    const int n = A.Size1;
    Q = std::move(A);
    R = TMatrix<double>(n, n);
    R.Nullify();
    for (int j = 0; j < n; ++j) {
        double* jthRColumn = R.GetColumn(j);
        double* jthQColumn = Q.GetColumn(j);
        double* innerProdVector = new double[blockSize];
        for (int i = 0; i < j; i += blockSize) {
            int actualBlockSize = std::min(blockSize, j - i);
            for (int k = 0; k < actualBlockSize; ++k) {
                innerProdVector[k] = 0;
            }
            double* ithQColumn = Q.GetColumn(i);
            for (int row = 0; row < n; ++row) {
                double coeff = jthQColumn[row];
                for (int k = 0; k < actualBlockSize; ++k) {
                    innerProdVector[k] += coeff * ithQColumn[k * n + row];
                }
            }

            for (int l = 0; l < actualBlockSize; ++l) {
                jthRColumn[i + l] = innerProdVector[l];
                for (int k = 0; k < n; ++k) {
                    jthQColumn[k] -= ithQColumn[l * n + k] * innerProdVector[l];
                }
            }
        }
        delete[] innerProdVector;
        double rjj = 0;
        for (int i = 0; i < n; ++i) {
            rjj += jthQColumn[i] * jthQColumn[i];
        }
        rjj = std::sqrt(rjj);
        for (int k = 0; k < n; ++k) {
            jthQColumn[k] /= rjj;
        }
        jthRColumn[j] = rjj;
    }
    return true;
}

bool QRDecompositionBlockOptimizedBlas(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R) {
    if (A.Size1 != A.Size2) {
        return false;
    }
    const int blockSize = 512;
    const int n = A.Size1;
    const int increment = 1;
    const double alpha = 1;
    const double beta = 0;
    const char trans = 't';
    Q = std::move(A);
    R = TMatrix<double>(n, n);
    R.Nullify();
    double* innerProdVector = new double[blockSize];
    for (int j = 0; j < n; ++j) {
        double* jthRColumn = R.GetColumn(j);
        double* jthQColumn = Q.GetColumn(j);
        for (int i = 0; i < j; i += blockSize) {
            int actualBlockSize = std::min<int>(blockSize, j - i);
            double* ithQColumn = Q.GetColumn(i);
            dgemv_(&trans, &n, &actualBlockSize, &alpha, ithQColumn, &n, jthQColumn, &increment, &beta, innerProdVector, &increment);

            for (int l = 0; l < actualBlockSize; ++l) {
                jthRColumn[i + l] = innerProdVector[l];
                innerProdVector[l] *= -1;
                daxpy_(&n, &innerProdVector[l], &ithQColumn[l * n], &increment, jthQColumn, &increment);
            }
        }
        double rjj = dnrm2_(&n, jthQColumn, &increment);
        jthRColumn[j] = rjj;
        const double coeff = 1 / rjj;
        dscal_(&n, &coeff, jthQColumn, &increment);
    }
    delete[] innerProdVector;
    return true;
}