/**
 *  A benchmark executable for the construction of the Hubbard Hamiltonian matrix
 */

#include <benchmark/benchmark.h>

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"


static void constructHamiltonian(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::Hubbard hubbard (fock_space);

    GQCP::HamiltonianParameters ham_par = GQCP::HamiltonianParameters::Random(K);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::MatrixXd hamiltonian = hubbard.constructHamiltonian(ham_par);

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
