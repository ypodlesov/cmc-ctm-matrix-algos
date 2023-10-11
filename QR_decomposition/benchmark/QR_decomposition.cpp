#include <benchmark/benchmark.h>
#include <QR_decomposition.h>

static void BM_QRDecomposition(benchmark::State& state) {
  for (auto _ : state) {
      for (auto _ : state) {
          state.PauseTiming();
          auto A = TMatrix<double>::CreateRandom(state.range(0), state.range(0));
          TMatrix<double> Q, R;
          state.ResumeTiming();
          QRDecomposition(A, Q, R);
      }
  }
}

BENCHMARK(BM_QRDecomposition)->Arg(1 << 8)->Arg(1 << 9)->Arg(1 << 10)->Arg(1 << 11)->Arg(1 << 12);

BENCHMARK_MAIN();