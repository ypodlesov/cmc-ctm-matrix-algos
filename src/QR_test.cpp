#include <gtest/gtest.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <ctime>

#include "QR.h"
namespace {
    using boost::random::mt19937;
    using boost::random::uniform_int_distribution; 
    const int lower_bound = -1e6;
    const int upper_bound = 1e6;
}

class QRTestBase : public ::testing::Test {
public:
    QRTestBase() 
        : N{1000}
        , A(N, N)
        , generator(std::time(0))
        , range(lower_bound, upper_bound)
    {
    }
    virtual void SetUp() override {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                A(i, j) = range(generator);
            }
        }
        B = A;
    } 

protected:
    double MatrixDist(const matrix<double>& A, const matrix<double>& B) {
        double delta = 0.0;
        for (size_t i = 0; i < A.size1(); ++i) {
            for (size_t j = 0; j < A.size2(); ++j) {
                delta += abs(A(i, j) - B(i, j));
            }
        }
        return delta;
    }

    const size_t N;
    matrix<double> A;
    matrix<double> B;
    matrix<double> Q;
    matrix<double> R;
    mt19937 generator;
    uniform_int_distribution<> range;
};

TEST_F(QRTestBase, QRTest) {
    ASSERT_TRUE(QR_decomposition(B, Q, R));
    matrix<double> QR = prod(Q, R);

    ASSERT_TRUE(QR.size1() == A.size1() && QR.size2() == A.size2());
    double delta = MatrixDist(A, QR);
    ASSERT_TRUE(rough_eq(delta, 0.0, 1e-3));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}