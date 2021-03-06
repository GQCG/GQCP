/**
 *  A benchmark executable for the construction of the Hubbard Hamiltonian matrix. The number of sites is kept at 8, while the number of electron pairs varies from 2 to 4.
 */

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 5; ++i) {  // Needs an `int` instead of a `size_t`.
        b->Args({8, i});           // The number of sites, the number of electron pairs.
    }
}


static void constructHamiltonian(benchmark::State& state) {

    // Prepare a random Hubbard model Hamiltonian.
    const size_t K = state.range(0);    // The number of lattice sites.
    const size_t N_P = state.range(1);  // The number of electron pairs.
    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


    // Set up the full spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto H = onv_basis.evaluateOperatorDense(hubbard_hamiltonian);

        benchmark::DoNotOptimize(H);  // Make sure that the variable is not optimized away by compiler.
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
