#include <gtest/gtest.h>
#include <math.h>
#include <qr_decomposition.h>
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
        A = TMatrix<double>::CreateRandom(N, N);
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

INSTANTIATE_TEST_SUITE_P(TQRTest, TQRTestBase, testing::Values(4, 64, 128, 256));