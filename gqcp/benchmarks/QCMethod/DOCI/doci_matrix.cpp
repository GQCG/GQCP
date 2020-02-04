/**
 *  A benchmark executable for the construction of the DOCI Hamiltonian
 */

#include <benchmark/benchmark.h>

#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


static void constructHamiltonian(benchmark::State& state) {

    // Prepare the Hamiltonian
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ONVBasis fock_space (K, N);
    GQCP::DOCI doci (fock_space);
    auto hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::SquareMatrix<double> H = doci.constructHamiltonian(hamiltonian);

        benchmark::DoNotOptimize(H);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {
        b->Args({16, i});  // orbitals, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
