/**
 *  A benchmark executable for the construction of the FCI Hamiltonian, which is constructed for full spin-resolved SpinUnresolvedONV bases for 8 orbitals and 2 to 6 electron pairs.
 */

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // Needs an `int` instead of a `size_t`.
        b->Args({8, i});           // The number of spatial orbitals, the number of electron pairs.
    }
}


static void constructHamiltonian(benchmark::State& state) {

    const size_t K = state.range(0);    // The number of spatial orbitals.
    const size_t N_P = state.range(1);  // The number of electron pairs.


    // Prepare the second-quantized Hamiltonian and set up the full spin-resolved ONV basis.
    // Note that the Hamiltonian is not necessarily expressed in an orthonormal basis, but this doesn't matter here.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


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
