/*
 *  A benchmark executable to check the performance of finding FCI results using a dense and iterative (Davidson) algorithm.
 *
 *  The system of interest is a linear H-chain composed of 4 to 10 hydrogen atoms.
 */

#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <benchmark/benchmark.h>


static void CustomArguments(benchmark::internal::Benchmark* b) {
    for (int i = 4; i < 11; i++) {  // Needs an `int` instead of a `size_t`.
        b->Args({i, 4});            // The number of hydrogen nuclei, 4 electrons.
    }
}


/**
 *  A dense benchmark.
 */
static void fci_dense_molecule(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const size_t N = state.range(1);  // The number of electrons.
    const auto N_P = N / 2;           // The number of electron pairs.
    const auto charge = static_cast<int>(number_of_H_atoms - N);


    // Set up the molecular Hamiltonian in the canonical RHF basis.
    // Construct the initial spinor basis.
    const auto molecule = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // Represented in AO basis.

    // Solve the SCF equations using a plain solver to find the canonical spinors.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::transform(rhf_parameters.expansion(), spinor_basis, hamiltonian);


    // Do the FCI calculation by setting up a full spin-resolved ONV basis, an eigenvalue problem solver and a corresponding environment.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    auto environment = GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();


    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

        benchmark::DoNotOptimize(electronic_energy);  // Make sure that the variable is not optimized away by the compiler.
    }

    state.counters["Hydrogen nuclei"] = number_of_H_atoms;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = onv_basis.dimension();
}


/**
 *  An iterative (Davidson) benchmark.
 */
static void fci_davidson_molecule(benchmark::State& state) {

    const auto number_of_H_atoms = state.range(0);
    const size_t N = state.range(1);  // The number of electrons.
    const auto N_P = N / 2;           // The number of electron pairs.
    const auto charge = static_cast<int>(number_of_H_atoms - N);


    // Set up the molecular Hamiltonian in the canonical RHF basis.
    // Construct the initial spinor basis.
    const auto molecule = GQCP::Molecule::HChain(number_of_H_atoms, 0.742, charge);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // Represented in AO basis.


    // Solve the SCF equations using a plain solver to find the canonical spinors.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::transform(rhf_parameters.expansion(), spinor_basis, hamiltonian);


    // Do the FCI calculation by setting up a full spin-resolved ONV basis, an eigenvalue problem solver and a corresponding environment.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    const auto initial_guess = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();


    // Code inside this loop is measured repeatedly.
    for (auto _ : state) {
        const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

        benchmark::DoNotOptimize(electronic_energy);  // Make sure that the variable is not optimized away by the compiler.
    }

    state.counters["Hydrogen nuclei"] = K;
    state.counters["Electrons"] = N;
    state.counters["Dimension"] = onv_basis.dimension();
}


BENCHMARK(fci_davidson_molecule)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK(fci_dense_molecule)->Unit(benchmark::kMillisecond)->Apply(CustomArguments);
BENCHMARK_MAIN();
