#include "matrix.h"
#include "vector.h"

extern "C" {
    double dnrm2_(const int *, const double *, const int *);
    double ddot_(const int *, const double *, const int *, const double *, const int *);
    void dscal_(const int *, const double *, const double *, const int *);
    void daxpy_(const int *, const double *, const double *, const int *, const double *, const int *);
    void dgemv_(const char &, const int &, const int &, const double &, const double *, const int &, const double *, const int &, const double &, double *, const int &);
    void dgemm_(const char &, const char &, 
        const int &, const int &, const int &, 
        const double &, const double *, const int &, const double *, const int &, 
        const double &, double *, const int &);
}

bool QRDecomposition(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);
bool QRDecompositionBlas(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);
bool QRDecompositionBlockOptimized(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);
bool QRDecompositionBlockOptimizedBlas(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);