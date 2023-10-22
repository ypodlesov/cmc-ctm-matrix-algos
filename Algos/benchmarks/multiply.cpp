#include <benchmark/benchmark.h>
#include <matrix.h>

static void BM_Prod(benchmark::State& state) {
    for (auto _ : state) {
        for (auto _ : state) {
            state.PauseTiming();
            auto A = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
            auto B = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
            state.ResumeTiming();
            TMatrix<double>::Prod(A, B);
        }
    }
}

BENCHMARK(BM_Prod)->Arg(1 << 6)->Arg(1 << 7)->Arg(1 << 8)->Arg(1 << 9)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

static void BM_BlockProd(benchmark::State& state) {
    for (auto _ : state) {
        for (auto _ : state) {
            state.PauseTiming();
            auto A = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
            auto B = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
            state.ResumeTiming();
            TMatrix<double>::BlockProd(A, B, {32, 32, 32});
        }
    }
}

BENCHMARK(BM_BlockProd)->Arg(1 << 6)->Arg(1 << 7)->Arg(1 << 8)->Arg(1 << 9)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

BENCHMARK_MAIN();