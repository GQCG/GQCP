/**
 *  A benchmark executable for DOCI calculations on CO
 */

#include <benchmark/benchmark.h>

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


static void test_case(benchmark::State& state) {
    // Benchmark a DOCI calculation based on a given FCIDUMP file

    // Prepare the Hamiltonian
    auto hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("../../benchmarks/benchmark_input/co_631g_klaas.FCIDUMP");

    // The species contains 14 electrons and 28 basis functions, this requires a single Fock Space of 28 orbitals and 7 electrons
    GQCP::FockSpace fock_space (hamiltonian.dimension(), 7);  // dim = 1184040

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::VectorX<double> initial_guess = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions solver_options (initial_guess);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (doci, hamiltonian);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Orbitals"] = hamiltonian.dimension();
    state.counters["Electron pairs"] = 7;
    state.counters["Dimension"] = fock_space.get_dimension();
}


// Perform the benchmarks
BENCHMARK(test_case)->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();
