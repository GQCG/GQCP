/**
 *  A benchmark executable for the diagonalization of Hubbard Hamiltonian matrices
 */

#include <benchmark/benchmark.h>

#include "CISolver/CISolver.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"


static void constructHamiltonian(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::Hubbard hubbard (fock_space);

    GQCP::HamiltonianParameters ham_par = GQCP::HamiltonianParameters::Random(K);
    GQCP::DenseSolverOptions solver_options;

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (hubbard, ham_par);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Sites"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 10; ++i) {  // need int instead of size_t
        b->Args({i, 2});  // sites, electron pairs
    }
}


// Perform the benchmarks
BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
