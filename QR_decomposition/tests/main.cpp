// #include "multiply_block.cpp"
// #include "multiply_parallel.cpp"
#include <qr_decomposition.h>
#include <chrono>

int main(int argc, char** argv) {
    std::srand(std::time(nullptr));
    const size_t N = 4096;
    auto A = TMatrix<double>::CreateRandom(N, N);
    auto B = A;
    TMatrix<double> Q, R;
    auto startTime = std::chrono::steady_clock::now();
    std::cout << QRDecompositionClassicBlas(A, Q, R) << std::endl;
    std::cout << std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime) << std::endl;
    // testing::InitGoogleTest(&argc, argv);
    // return RUN_ALL_TESTS();
}