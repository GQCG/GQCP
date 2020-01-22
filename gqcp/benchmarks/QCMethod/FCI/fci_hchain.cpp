/**
 *  A benchmark executable for the FCI
 */

#include <benchmark/benchmark.h>

#include "Basis/transform.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrix.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"



/**
 *  DAVIDSON
 */
static void fci_davidson_hchain(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const auto number_of_electrons = state.range(1);
    const auto charge = static_cast<int>(number_of_H_atoms - number_of_electrons);


    // Create the molecular Hamiltonian for this molecule and basis
    const auto hchain = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in AO basis


    // Do an RHF SCF procedure
    const auto K = sq_hamiltonian.dimension();
    const auto N_P = hchain.numberOfElectrons()/2;

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(hchain.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrix<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Diagonalize the FCI Hamiltonian in the RHF basis
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());
    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    GQCP::VectorX<double> initial_guess = fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions solver_options (initial_guess);


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (fci, sq_hamiltonian);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen nuclei"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}


/**
 *  DENSE
 */
static void fci_dense_hchain(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const auto number_of_electrons = state.range(1);
    const auto charge = static_cast<int>(number_of_H_atoms - number_of_electrons);


    // Create the molecular Hamiltonian for this molecule and basis
    const auto hchain = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);
    auto K = sq_hamiltonian.dimension();
    auto N_P = hchain.numberOfElectrons()/2;

    // Create a plain RHF SCF solver and solve the SCF equations
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(hchain.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrix<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());
    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    GQCP::DenseSolverOptions solver_options;

    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (fci, sq_hamiltonian);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen nuclei"] = K;
    state.counters["Electrons"] = 2*N_P;
    state.counters["Dimension"] = fock_space.get_dimension();
}


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 11; i++) {  // need int instead of size_t
        b->Args({i, 4});  // number of hydrogen nuclei, 4 electrons
    }
}


// Perform the benchmarks
BENCHMARK(fci_davidson_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK(fci_dense_hchain)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
