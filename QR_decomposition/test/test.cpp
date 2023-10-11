#include <cstdlib>
#include <ctime>
#include <gtest/gtest.h>
#include <math.h>
#include <QR_decomposition.h>
#include <common.h>
#include <iostream>

class TQRTestBase : public testing::TestWithParam<size_t> {
public:
    TQRTestBase() 
        : N{GetParam()}
    {
    }

    virtual void SetUp() override {
        std::srand(std::time(nullptr));
        A = TMatrix<double>(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                A(i, j) = std::rand() % N;
            }
        }
        B = A;
    } 

protected:
    const size_t N;
    TMatrix<double> A;
    TMatrix<double> B;
    TMatrix<double> Q;
    TMatrix<double> R;
};

TEST_P(TQRTestBase, QRDecompositionTest) {
    ASSERT_TRUE(QRDecomposition(B, Q, R));
    TMatrix C(TMatrix<double>::Prod(Q, R));
    ASSERT_EQ(A, C);
    ASSERT_TRUE(TMatrix<double>::IsTriangular(R, ETriangularType::Upper));
    ASSERT_EQ(TMatrix<double>::Prod(Q, Q.Transpose()), TMatrix<double>::CreateIdentityMatrix(N));
}

INSTANTIATE_TEST_SUITE_P(TQRTestInstance, TQRTestBase, testing::Values(3, 64, 128, 256));

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}