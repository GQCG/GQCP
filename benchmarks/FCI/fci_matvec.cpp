/**
 *  A benchmark executable for the FCI matvec
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/FCI.hpp"


static void matvec(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::FCI fci (fock_space);

    // Random ham_par
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    Eigen::VectorXd diagonal = fci.calculateDiagonal(ham_par);
    Eigen::VectorXd random = Eigen::VectorXd::Random(diagonal.rows());

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::VectorXd matvec = fci.matrixVectorProduct(ham_par, random, diagonal);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(matvec);
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i){
        b->Args({10,i});
    }
}

// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
