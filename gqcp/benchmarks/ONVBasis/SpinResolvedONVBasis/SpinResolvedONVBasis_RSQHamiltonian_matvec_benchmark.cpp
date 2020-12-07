/**
 *  A benchmark executable that times the performance of one FCI matrix-vector product in a full spin-resolved ONV basis with 10 orbitals and 2 to 5 electron pairs.
 */

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({10, i});          // spatial orbitals, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    const size_t K = state.range(0);    // number of spatial orbitals
    const size_t N_P = state.range(1);  // number of electron pairs


    // Prepare the second-quantized Hamiltonian and set up the full spin-resolved ONV basis.
    // Note that the Hamiltonian is not necessarily expressed in an orthonormal basis, but this doesn't matter here.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    const auto x = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis).coefficients();

    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto matvec = onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, x);

        benchmark::DoNotOptimize(matvec);  // Make sure that the variable is not optimized away by compiler.
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
