/**
 *  A benchmark executable for DOCI calculations on CO in a 6-31G basisset.
 */

#include <benchmark/benchmark.h>

#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


static void test_case(benchmark::State& state) {

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/co_631g_klaas.FCIDUMP");

    // Do the DOCI calculation by setting up a full seniority-zero ONV basis, a Davidson solver and a corresponding environment.
    // The system contains 14 electrons and 28 basis functions.
    const auto K = sq_hamiltonian.dimension();
    const auto N_P = 7;

    const GQCP::SeniorityZeroONVBasis onv_basis (K, N_P);  // dim = 1184040
    const GQCP::VectorX<double> initial_guess = onv_basis.hartreeFockExpansion();

    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Davidson();


    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto electronic_energy = GQCP::QCMethod::CI(onv_basis).optimize(solver, environment).groundStateEnergy();

        benchmark::DoNotOptimize(electronic_energy);  // make sure that the variable is not optimized away by compiler
    }

    state.counters["Spatial orbitals"] = K;
    state.counters["Electron pairs"] = N_P;
    state.counters["Dimension"] = onv_basis.get_dimension();  // dim = 1184040
}


BENCHMARK(test_case)->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();
