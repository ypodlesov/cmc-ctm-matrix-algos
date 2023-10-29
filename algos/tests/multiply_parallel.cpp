#include <gtest/gtest.h>
#include <matrix.h>
#include <math.h>
#include <common.h>
#include <iostream>

TEST_P(TMultiplyTestBase, ParallelMultiplyTest) {
    TMatrix C(TMatrix<double>::Prod(A, B));
    TMatrix D(TMatrix<double>::ParallelProd(A, B));
    ASSERT_EQ(C, D);
}

INSTANTIATE_TEST_SUITE_P(TParallelMultiplyTest, TMultiplyTestBase, testing::Values(4, 64, 128, 256));