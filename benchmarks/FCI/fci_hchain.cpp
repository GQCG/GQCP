/**
 *  A benchmark executable for the FCI
 */

#include <benchmark/benchmark.h>

#include "CISolver/CISolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Molecule.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"



/**
 *  DAVIDSON
 */
static void fci_davidson_hchain(benchmark::State& state) {

    int64_t number_of_H_atoms = state.range(0);
    int64_t number_of_electrons = state.range(1);
    auto charge = static_cast<int>(number_of_H_atoms - number_of_electrons);

    GQCP::Molecule hchain = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(hchain, "STO-3G");
    auto K = mol_ham_par.get_K();
    auto N_P = hchain.get_N()/2;
    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, hchain);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    mol_ham_par.transform(rhf.get_C());
    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd initial_guess = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_guess);

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (fci, mol_ham_par);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen atoms"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}

/**
 *  DENSE
 */
static void fci_dense_hchain(benchmark::State& state) {

    int64_t number_of_H_atoms = state.range(0);
    int64_t number_of_electrons = state.range(1);
    auto charge = static_cast<int>(number_of_H_atoms - number_of_electrons);

    GQCP::Molecule hchain = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(hchain, "STO-3G");
    auto K = mol_ham_par.get_K();
    auto N_P = hchain.get_N()/2;
    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, hchain);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    mol_ham_par.transform(rhf.get_C());
    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    numopt::eigenproblem::DenseSolverOptions solver_options;

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (fci, mol_ham_par);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen atoms"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 11; i++) {  // need int instead of size_t
        b->Args({i, 4});  // number of hydrogen atoms, 4 electrons
    }
}


// Perform the benchmarks
BENCHMARK(fci_davidson_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK(fci_dense_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
