#include <gtest/gtest.h>
#include <matrix.h>
#include <math.h>
#include <common.h>
#include <iostream>
#include "multiply_base.h"

TEST_P(TMultiplyTestBase, BlockMultiplyTest) {
    TMatrix C(TMatrix<double>::Prod(A, B));
    TMatrix D(TMatrix<double>::BlockProd(A, B, {32, 32, 32}));
    ASSERT_EQ(C, D);
}

INSTANTIATE_TEST_SUITE_P(TBlockMultiplyTest, TMultiplyTestBase, testing::Values(4, 64, 128, 256));