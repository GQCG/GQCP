/**
 *  A benchmark executable for the Hubbard constructHamiltonian
 */

#include "benchmark/benchmark.h"

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"


static void constructHamiltonian(benchmark::State& state) {
    // Prepare parameters
    size_t K = state.range(0);
    size_t N = state.range(1);
    GQCP::ProductFockSpace fock_space (K, N, N);
    GQCP::Hubbard hubbard (fock_space);

    // Random ham_par
    GQCP::HamiltonianParameters ham_par = GQCP::constructRandomHamiltonianParameters(K);
    numopt::eigenproblem::DenseSolverOptions solver_options;

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (hubbard, ham_par);
        ci_solver.solve(solver_options);
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(ci_solver);
    }

    state.counters["Orbitals"] = K;
    state.counters["Electron pairs"] = N;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 2; i < 4; ++i){
          b->Args({6,i});
    }
}

// Perform the benchmarks
BENCHMARK(constructHamiltonian)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
