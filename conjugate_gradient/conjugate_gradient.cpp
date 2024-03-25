#include "conjugate_gradient.h"
#include <helpers.h>

bool ConjugateGradient(const TMatrix<double>& A, const TVector<double>& b, TVector<double>& x) {
    if (!A || !b || A.Size2 != b.Size) {
        return false;
    }
    TVector<double> currentResidual, nextResidual, currentP, nextP, currentX, nextX;
    currentResidual = currentP = b;
    currentX.Nullify();
    nextX.Nullify();
    double currentAlpha, currentBeta;
    for (size_t j = 0; j < 100 && !RoughEq<double>(TVector<double>::Norm2(currentResidual), 0); ++j) {
        TVector<double> Ap;
        TMatrix<double>::MVMultiply(A, currentP, Ap);
        double innerProduct = TVector<double>::InnerProd(Ap, currentP);
        if (RoughEq(innerProduct, 0)) {
            x = std::move(currentX);
            return true;
        }
        double currentResidualNorm = TVector<double>::InnerProd(currentResidual, currentResidual);
        currentAlpha = currentResidualNorm / innerProduct;
        nextX = currentX + currentAlpha * currentP;
        nextResidual = currentResidual - currentAlpha * Ap;
        if (RoughEq(currentResidualNorm, 0)) {
            x = std::move(currentX);
            return true; 
        }
        currentBeta = TVector<double>::InnerProd(nextResidual, nextResidual) / currentResidualNorm;
        nextP = nextResidual + currentBeta * currentP;
    }
    return true;
}