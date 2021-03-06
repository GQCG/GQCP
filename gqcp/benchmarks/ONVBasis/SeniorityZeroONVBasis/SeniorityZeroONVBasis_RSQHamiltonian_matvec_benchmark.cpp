/**
 *  A benchmark executable for the DOCI matrix-vector product. The system of interest has K=28 spatial orbitals and N_P=5-8 electron pairs.
 */

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {  // need int instead of size_t
        b->Args({28, i});          // spatial orbitals, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    // Set up a random restricted SQHamiltonian and a seniority-zero ONV basis.
    const size_t K = state.range(0);    // The number of spatial orbitals.
    const size_t N_P = state.range(1);  // The number of electron pairs.

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);  // This Hamiltonian is not necessarily expressed in an orthonormal basis, but this doesn't matter here.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, N_P};

    const auto x = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::Random(onv_basis).coefficients();

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
