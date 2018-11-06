/**
 *  A benchmark executable for the DOCI matvec
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/DOCI.hpp"


static void matvec(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::FockSpace fock_space (K, N);
    GQCP::DOCI doci (fock_space);

    // Random ham_par
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(ham_par);
    Eigen::VectorXd random = Eigen::VectorXd::Random(diagonal.rows());

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::VectorXd matvec = doci.matrixVectorProduct(ham_par, random, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }

    state.counters["Electron pairs"] = N;
    state.counters["Orbitals"] = K;
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 10; i <= 20; ++i){
          b->Args({i, i/2});
    }
}

// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();