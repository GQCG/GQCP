/**
 *  A benchmark executable for the one electron matvec performance of FockSpace
 */

#include <benchmark/benchmark.h>

#include "FockSpace/FockSpace.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"


static void matvec(benchmark::State& state) {

    // Prepare the Fock space
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);

    GQCP::QCMatrix<double> one_op_par = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::ScalarSQOneElectronOperator<double> sq_one_op {one_op_par};
    GQCP::VectorX<double> diagonal = fock_space.evaluateOperatorDiagonal(sq_one_op);
    GQCP::VectorX<double> x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::VectorX<double> matvec = fock_space.evaluateOperatorMatrixVectorProduct(sq_one_op, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 10; ++i) {  // need int instead of size_t
        b->Args({28, i});  // orbitals, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
