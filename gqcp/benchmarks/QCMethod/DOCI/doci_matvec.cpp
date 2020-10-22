/**
 *  A benchmark executable for the DOCI matrix-vector product. The system of interest has K=28 spatial orbitals and N_P=5-8 electron pairs.
 */

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {  // need int instead of size_t
        b->Args({28, i});          // spatial orbitals, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    // Set up a random second-quantized Hamiltonian and a doubly-occupied ONV basis.
    const size_t K = state.range(0);
    const size_t N_P = state.range(1);  // number of electron pairs

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Random(K);  // not necessarily in a non-orthonormal basis, but this doesn't matter here
    const GQCP::SeniorityZeroONVBasis onv_basis {K, N_P};

    GQCP::DOCI doci_builder {onv_basis};  // the DOCI HamiltonianBuilder

    const auto diagonal = doci_builder.calculateDiagonal(sq_hamiltonian);
    const auto x = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::Random(onv_basis).coefficients();

    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto matvec = doci_builder.matrixVectorProduct(sq_hamiltonian, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
