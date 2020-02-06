/**
 *  A benchmark executable for the DOCI matvec
 */

#include <benchmark/benchmark.h>

#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


/**
 *  Benchmark the performance of a matrix-vector multiplication for DOCI, using K=28 spatial orbitals and N_P=5-8 electron pairs.
 */
static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 5; i < 9; ++i) {  // need int instead of size_t
        b->Args({28, i});  // spatial orbitals, electron pairs
    }
}


static void matvec(benchmark::State& state) {

    // Set up a random SQHamiltonian and a doubly occupied ONV basis
    const size_t K = state.range(0);
    const size_t N_P = state.range(1);

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);
    GQCP::ONVBasis doubly_occupied_onv_basis (K, N_P);

    GQCP::DOCI doci (doubly_occupied_onv_basis);

    const auto diagonal = doci.calculateDiagonal(sq_hamiltonian);
    const auto x = doubly_occupied_onv_basis.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto matvec = doci.matrixVectorProduct(sq_hamiltonian, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = doubly_occupied_onv_basis.get_dimension();
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
