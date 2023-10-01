#include <iostream>
#include "QR.h"

int main() {

    matrix<double> A(3, 3);
    for (unsigned i = 0; i < A.size1 (); ++ i) {
        for (unsigned j = 0; j < A.size2 (); ++ j) {
            A (i, j) = 3 * i + j;
        }
    }
    std::cout << A << std::endl;
    auto B = A;
    matrix<double> Q;
    triangular_matrix<double, upper> R;
    if (QR_decomposition(B, Q, R)) {
        std::cout << prod(Q, R) << std::endl;
    } else {
        std::cout << "unable to process QR decomposition" << std::endl;
    }

}