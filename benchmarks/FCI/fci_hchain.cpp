/**
 *  A benchmark executable for the FCI
 */

#include "benchmark/benchmark.h"

#include "CISolver/CISolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "Molecule.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

/**
 *  AUXILIARY FUNCTIONS
 */

/**
 *  @param  chain_length
 *  @param  interatomic_distance
 *  @param  electrons
 *
 *  @return hydrogen-chain with specified length, spacing and electrons
 */
static GQCP::Molecule makeHChain(size_t chain_length, double interatomic_distance, size_t electrons){
    if (chain_length == 0) {
        throw std::invalid_argument("Can not return a molecule consisting of zero atoms");
    }

    std::vector<GQCP::Atom> hs;
    size_t half_num = chain_length/2;
    double outer_x_pos = -(half_num*interatomic_distance);
    if (!(chain_length & 2)) {
        outer_x_pos -= interatomic_distance/2;
    }

    for (size_t i = 0; i < chain_length; i++) {
        double x_pos = outer_x_pos + i*interatomic_distance;
        hs.emplace_back(1, x_pos, 0, 0);
    }

    int charge = chain_length-electrons;
    return GQCP::Molecule(hs, charge);
}

/**
 *  BENCHMARKS
 */

/**
 *  DAVIDSON
 */
static void fci_davidson_hchain(benchmark::State& state) {
    GQCP::Molecule hchain = makeHChain(state.range(0), 0.742 , state.range(1));
    auto ao_basis = std::make_shared<GQCP::AOBasis>(hchain, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
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
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(ci_solver);
    }

    state.counters["Hydrogen atoms"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}

/**
 *  DENSE
 */
static void fci_dense_hchain(benchmark::State& state) {
    GQCP::Molecule hchain = makeHChain(state.range(0), 0.742 , state.range(1));
    auto ao_basis = std::make_shared<GQCP::AOBasis>(hchain, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
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
        // Make sure the variable is not optimized away by compiler
        benchmark::DoNotOptimize(ci_solver);
    }

    state.counters["Hydrogen atoms"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}

static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 11; i++) { // Hydrogen atoms
        // 4 electrons

        b->Args({i, 4});

    }
}

// Perform the benchmarks
BENCHMARK(fci_davidson_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK(fci_dense_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);


BENCHMARK_MAIN();
