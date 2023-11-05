#include "matrix.h"
#include "vector.h"

bool QRDecomposition(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);
bool QRDecompositionBlockOptimized(TMatrix<double>& A, TMatrix<double>& Q, TMatrix<double>& R);