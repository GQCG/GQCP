/**
 *  A benchmark executable for the DOCI matvec
 */

#include <benchmark/benchmark.h>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "FockSpace/FockSpace.hpp"


static void matvec(benchmark::State& state) {

    // Prepare the Hamiltonian
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);

    GQCP::QCMatrix<double> one_electron_parameters = GQCP::QCMatrix<double>::Random(K, K);
    GQCP::ScalarSQOneElectronOperator<double> sq_one_electron ({one_electron_parameters});
    GQCP::VectorX<double> diagonal = fock_space.evaluateOperatorDiagonal(sq_one_electron);
    GQCP::VectorX<double> x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::VectorX<double> matvec = fock_space.evaluateOperatorMatvec(sq_one_electron, x, diagonal);

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
