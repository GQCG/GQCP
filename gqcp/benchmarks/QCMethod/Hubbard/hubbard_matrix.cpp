/**
 *  A benchmark executable for the construction of the Hubbard Hamiltonian matrix
 */

#include <benchmark/benchmark.h>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"


static void constructHamiltonian(benchmark::State& state) {

    // Prepare the Hamiltonian
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::Hubbard hubbard (fock_space);

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::SquareMatrix<double> hamiltonian = hubbard.constructHamiltonian(sq_hamiltonian);

        benchmark::DoNotOptimize(hamiltonian);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 5; ++i) {  // need int instead of size_t
        b->Args({8, i});  // sites, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
