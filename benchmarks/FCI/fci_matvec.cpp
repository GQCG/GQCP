/**
 *  A benchmark executable for the FCI matvec
 */

#include <benchmark/benchmark.h>

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/FCI.hpp"


static void matvec(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::FCI fci (fock_space);

    GQCP::HamiltonianParameters<double> ham_par = GQCP::HamiltonianParameters<double>::Random(K);
    Eigen::VectorXd diagonal = fci.calculateDiagonal(ham_par);
    Eigen::VectorXd x = fock_space.randomExpansion();

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        Eigen::VectorXd matvec = fci.matrixVectorProduct(ham_par, x, diagonal);

        benchmark::DoNotOptimize(matvec);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 6; ++i) {  // need int instead of size_t
        b->Args({10, i});  // orbitals, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(matvec)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
