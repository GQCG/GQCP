/**
 *  A benchmark executable that tests the one-electron matrix-vector product for a full spin-unresolved ONV basis. The number of spinors is 28, the number of electrons varies from 5 to 9.
 */

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 10; ++i) {  // Needs an `int` instead of a `size_t`.
        b->Args({28, i});           // The number of spinors, the number of electrons.
    }
}


static void matvec(benchmark::State& state) {

    // Prepare a full spin-unresolved ONV basis and a random one-electron operator.
    const size_t M = state.range(0);  // The number of spinors.
    const size_t N = state.range(1);  // The number of electrons.
    GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    const auto f = GQCP::ScalarGSQOneElectronOperator<double>::Random(M);
    const auto x = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis).coefficients();

    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        GQCP::VectorX<double> matvec = onv_basis.evaluateOperatorMatrixVectorProduct(f, x);

        benchmark::DoNotOptimize(matvec);  // Make sure that the variable is not optimized away by compiler.
    }

    state.counters["Spinors"] = M;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
