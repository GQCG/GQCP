/**
 *  A benchmark for the construction of a random DOCI Hamiltonian matrix in a doubly-occupied spin-unresolved basis of 16 spatial orbitals and 5-8 electron pairs.
 */

#include <benchmark/benchmark.h>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {
        b->Args({16, i});  // spatial orbitals, electron pairs
    }
}


static void constructHamiltonian(benchmark::State& state) {

    // Set up a random SQHamiltonian and a doubly occupied SpinUnresolvedONV basis
    const size_t K = state.range(0);  // number of spatial orbitals
    const size_t N_P = state.range(1);  // number of electron pairs

    const auto hamiltonian = GQCP::SQHamiltonian<double>::Random(K);
    GQCP::SpinUnresolvedONVBasis doubly_occupied_onv_basis (K, N_P);

    GQCP::DOCI doci (doubly_occupied_onv_basis);


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto H = doci.constructHamiltonian(hamiltonian);

        benchmark::DoNotOptimize(H);  // make sure the variable is not optimized away by compiler
    }


    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = doubly_occupied_onv_basis.get_dimension();
}


BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
