/**
 *  A benchmark executable for the Hubbard matrix-vector product. The number of sites is kept at 12, while the number of electron pairs varies from 2 to 6.
 */

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 7; ++i) {  // need int instead of size_t
        b->Args({12, i});          // sites, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    // Prepare a random Hubbard model Hamiltonian.
    const size_t K = state.range(0);    // The number of lattice sites.
    const size_t N_P = state.range(1);  // The number of electron pairs.
    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


    // Set up the full spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    const auto x = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis).coefficients();


    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto matvec = onv_basis.evaluateOperatorMatrixVectorProduct(hubbard_hamiltonian, x);

        benchmark::DoNotOptimize(matvec);  // Make sure that the variable is not optimized away by compiler.
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
