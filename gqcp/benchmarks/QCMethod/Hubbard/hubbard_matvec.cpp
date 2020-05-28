/**
 *  A benchmark executable for the Hubbard matrix-vector product. The number of sites is kept at 12, while the number of electron pairs varies from 2 to 6.
 */

#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 7; ++i) {  // need int instead of size_t
        b->Args({12, i});          // sites, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    // Prepare a random Hubbard model Hamiltonian.
    const size_t K = state.range(0);    // number of sites
    const size_t N_P = state.range(1);  // number of electron pairs
    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


    // Set up the full spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    const GQCP::Hubbard hubbard_builder {onv_basis};
    const auto diagonal = hubbard_builder.calculateDiagonal(hubbard_hamiltonian);
    const auto x = onv_basis.randomExpansion();


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto matvec = hubbard_builder.matrixVectorProduct(hubbard_hamiltonian, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
