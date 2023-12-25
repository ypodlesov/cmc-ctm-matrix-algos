#include <benchmark/benchmark.h>
#include <qr_decomposition.h>

// static void BM_QRDecomposition(benchmark::State& state) {
//     for (auto _ : state) {
//         for (auto _ : state) {
//             state.PauseTiming();
//             auto A = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
//             TMatrix<double> Q, R;
//             state.ResumeTiming();
//             QRDecompositionBlas(A, Q, R);
//         }
//     }
// }

// BENCHMARK(BM_QRDecomposition)->Arg(1 << 6)->Arg(1 << 7)->Arg(1 << 8)->Arg(1 << 9)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

static void BM_QRDecompositionBlockOptimized(benchmark::State& state) {
    for (auto _ : state) {
        for (auto _ : state) {
            state.PauseTiming();
            auto A = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
            TMatrix<double> Q, R;
            state.ResumeTiming();
            QRDecompositionBlockOptimizedBlas(A, Q, R);
        }
    }
}

BENCHMARK(BM_QRDecompositionBlockOptimized)->Arg(1 << 6)->Arg(1 << 7)->Arg(1 << 8)->Arg(1 << 9)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

BENCHMARK_MAIN();