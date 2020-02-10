// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "FCI"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"


/**
 *  Check if we can reproduce the FCI energy for H2//6-31G**, with a dense solver. The reference is taken from Cristina (cfr. Ayers' lab).
 */
BOOST_AUTO_TEST_CASE ( FCI_H2_dense ) {

    const double reference_energy = -1.1651486697;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "6-31G**");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);  // dimension = 100


    // Create a dense solver and corresponding environment and put them together in the QCMethod
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if a random rotation has no effect on the sum of the diagonal elements of the FCI Hamiltonian matrix.
 */
BOOST_AUTO_TEST_CASE ( FCI_rotated_diagonal_sum ) {

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);  // dimension = 100


    // Evaluate the diagonal of the FCI Hamiltonian matrix.
    const auto diagonal = onv_basis.evaluateOperatorDiagonal(sq_hamiltonian);


    // Rotate the Hamiltonian with a random unitary matrix, evaluate the diagonal of the FCI Hamiltonian matrix in that basis, and check the result.
    sq_hamiltonian.randomRotate();
    const auto diagonal_rotated = onv_basis.evaluateOperatorDiagonal(sq_hamiltonian);
    BOOST_CHECK(std::abs(diagonal.sum() - diagonal_rotated.sum()) < 1.0e-10);
}


/**
 *  Check if we can reproduce the FCI energy for H2O//STO-3G**, with a dense solver. The reference is taken from Psi4 and GAMESS-US.
 */
BOOST_AUTO_TEST_CASE ( FCI_H2O_dense ) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);  // dimension = 100


    // Create a dense solver and corresponding environment and put them together in the QCMethod
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}
