#include "conjugate_gradient.h"
// #include "helpers.h"
#include "matrix.h"
// #include "vector.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <cstdint>

using Catch::Generators::random;
using Catch::Generators::take;

static void Test(const uint32_t /* n */) {
    auto matrix = Matrix<double>(10, 10);
    return;
}

TEST_CASE("Size 10") {
    Test(10);
}