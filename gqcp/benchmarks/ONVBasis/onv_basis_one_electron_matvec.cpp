/**
 *  A benchmark executable for the one electron matvec performance of ONVBasis
 */

#include <benchmark/benchmark.h>

#include "ONVBasis/ONVBasis.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"


static void matvec(benchmark::State& state) {

    // Prepare a full ONV basis
    const size_t M = state.range(0);  // the number of spinors
    const size_t N = state.range(1);  // the number of electrons
    GQCP::ONVBasis full_onv_basis (M, N);

    GQCP::QCMatrix<double> one_op_par = GQCP::QCMatrix<double>::Random(M, M);
    GQCP::ScalarSQOneElectronOperator<double> sq_one_op {one_op_par};
    GQCP::VectorX<double> diagonal = full_onv_basis.evaluateOperatorDiagonal(sq_one_op);
    GQCP::VectorX<double> x = full_onv_basis.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::VectorX<double> matvec = full_onv_basis.evaluateOperatorMatrixVectorProduct(sq_one_op, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Spinors"] = M;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = full_onv_basis.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 10; ++i) {  // need int instead of size_t
        b->Args({28, i});  // spinors, electrons
    }
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
