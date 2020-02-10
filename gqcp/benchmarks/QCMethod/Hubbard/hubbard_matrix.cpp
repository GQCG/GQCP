/**
 *  A benchmark executable for the construction of the Hubbard Hamiltonian matrix. The number of sites is kept at 8, while the number of electron pairs varies from 2 to 4.
 */

#include <benchmark/benchmark.h>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 5; ++i) {  // need int instead of size_t
        b->Args({8, i});  // sites, electron pairs
    }
}


static void constructHamiltonian(benchmark::State& state) {

    // Prepare a random Hamiltonian, but use the Hubbard matrix construction method. TODO: this has to be changed in a model Hamiltonian refactor
    const auto K = state.range(0);  // number of sites
    const auto N_P = state.range(1);  // number of electron pairs
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);


    // Set up the full spin-resolved ONV basis
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);
    GQCP::Hubbard hubbard (onv_basis);


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::SquareMatrix<double> hamiltonian = hubbard.constructHamiltonian(sq_hamiltonian);

        benchmark::DoNotOptimize(hamiltonian);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.get_dimension();
}


BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
