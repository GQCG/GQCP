/**
 *  A benchmark for the construction of a random DOCI Hamiltonian matrix in a doubly-occupied spin-unresolved basis of 16 spatial orbitals and 5-8 electron pairs.
 */

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {
        b->Args({16, i});  // spatial orbitals, electron pairs
    }
}


static void constructHamiltonian(benchmark::State& state) {

    // Set up a random restricted SQHamiltonian and a seniority-zero ONV basis.
    const size_t K = state.range(0);    // The number of spatial orbitals.
    const size_t N_P = state.range(1);  // The number of electron pairs.

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);  // This Hamiltonian is not necessarily expressed in an orthonormal basis, but this doesn't matter here.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, N_P};

    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto H = onv_basis.evaluateOperatorDense(hamiltonian);

        benchmark::DoNotOptimize(H);  // Make sure that the variable is not optimized away by compiler.
    }


    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
