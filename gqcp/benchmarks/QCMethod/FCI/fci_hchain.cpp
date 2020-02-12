/**
 *  A benchmark executable to check the performance of finding FCI results using a dense and iterative (Davidson) algorithm.
 * 
 *  The system of interest is a linear H-chain composed of 4 to 10 hydrogen atoms.
 */

#include <benchmark/benchmark.h>

#include "Basis/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"



static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 11; i++) {  // need int instead of size_t
        b->Args({i, 4});  // number of hydrogen nuclei, 4 electrons
    }
}



/**
 *  DENSE
 */
static void fci_dense_molecule(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const auto N = state.range(1);  // number of electrons
    const auto N_P = N / 2;  // number of electron pairs
    const auto charge = static_cast<int>(number_of_H_atoms - N);


    // Set up the molecular Hamiltonian in the canonical RHF basis.
    // Construct the initial spinor basis.
    const auto molecule = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in AO basis

    // Solve the SCF equations using a plain solver to find the canonical spinors.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Do the FCI calculation by setting up a full spin-resolved ONV basis, an eigenvalue problem solver and a corresponding environment.
    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        const auto electronic_energy = GQCP::QCMethod::CI(onv_basis).optimize(solver, environment).groundStateEnergy();

        benchmark::DoNotOptimize(electronic_energy);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen nuclei"] = number_of_H_atoms;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = onv_basis.get_dimension();
}



/**
 *  DAVIDSON
 */
static void fci_davidson_molecule(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const auto N = state.range(1);  // number of electrons
    const auto N_P = N / 2;  // number of electron pairs
    const auto charge = static_cast<int>(number_of_H_atoms - N);


    // Set up the molecular Hamiltonian in the canonical RHF basis.
    // Construct the initial spinor basis.
    const auto molecule = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in AO basis


    // Solve the SCF equations using a plain solver to find the canonical spinors.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());



    // Diagonalize the FCI Hamiltonian in the canonical RHF basis
    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);
    GQCP::FCI fci (onv_basis);

    GQCP::VectorX<double> initial_guess = onv_basis.hartreeFockExpansion();
    GQCP::DavidsonSolverOptions solver_options (initial_guess);


    // Code inside this loop is measured repeatedly
    for (auto _ : state) {
        GQCP::CISolver ci_solver (fci, sq_hamiltonian);
        ci_solver.solve(solver_options);

        benchmark::DoNotOptimize(ci_solver);  // make sure the variable is not optimized away by compiler
    }

    state.counters["Hydrogen nuclei"] = K;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = onv_basis.get_dimension();
}


BENCHMARK(fci_davidson_molecule)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK(fci_dense_molecule)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
