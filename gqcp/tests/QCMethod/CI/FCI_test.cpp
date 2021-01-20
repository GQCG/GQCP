// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "FCI"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check if we can reproduce the FCI energy for H2//6-31G**, with a dense solver. The reference is taken from Cristina (cfr. Ayers' lab).
 */
BOOST_AUTO_TEST_CASE(FCI_H2_dense) {

    const double reference_energy = -1.1651486697;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G**"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if a random rotation has no effect on the sum of the diagonal elements of the FCI Hamiltonian matrix.
 */
BOOST_AUTO_TEST_CASE(FCI_rotated_diagonal_sum) {

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.


    // Evaluate the diagonal of the FCI Hamiltonian matrix.
    const auto diagonal = onv_basis.evaluateOperatorDiagonal(sq_hamiltonian);


    // Rotate the Hamiltonian with a random unitary transformation, evaluate the diagonal of the FCI Hamiltonian matrix in that basis, and check the result.
    sq_hamiltonian.rotate(GQCP::RTransformation<double>::RandomUnitary(K));
    const auto diagonal_rotated = onv_basis.evaluateOperatorDiagonal(sq_hamiltonian);
    BOOST_CHECK(std::abs(diagonal.sum() - diagonal_rotated.sum()) < 1.0e-10);
}


/**
 *  Check if we can reproduce the FCI energy for H2O//STO-3G**, with a dense solver. The reference is taken from Psi4 and GAMESS-US.
 */
BOOST_AUTO_TEST_CASE(FCI_H2O_dense) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if we can reproduce the FCI energy for H2//6-31G**, with a Davidson solver. The reference is taken from Cristina (cfr. Ayers' lab).
 */
BOOST_AUTO_TEST_CASE(FCI_H2_Davidson) {

    const double reference_energy = -1.1651486697;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "6-31G**"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.

    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto x0 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::HartreeFock(onv_basis).coefficients();  // Supply an initial guess.
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, x0);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if we can reproduce the FCI energy for H2O//STO-3G**, with a Davidson solver. The reference is taken from Psi4 and GAMESS-US.
 */
BOOST_AUTO_TEST_CASE(FCI_H2O_Davidson) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis(K, N_P, N_P);  // The dimension of this ONV basis is 100.


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto x0 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::HartreeFock(onv_basis).coefficients();  // Supply an initial guess.
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, x0);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check that, for a test system H6//STO-3G, the dense FCI energy equals the Davidson FCI energy.
 */

BOOST_AUTO_TEST_CASE(FCI_H6_dense_vs_Davidson) {

    // Create the molecular Hamiltonian in an AO basis.
    const GQCP::Molecule molecule = GQCP::Molecule::HChain(6, 1.1);
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    sq_hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto dense_environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto dense_solver = GQCP::EigenproblemSolver::Dense();
    const auto dense_electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(dense_solver, dense_environment).groundStateEnergy();


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto x0 = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::HartreeFock(onv_basis).coefficients();  // Supply initial guess.
    auto davidson_environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, x0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson();
    const auto davidson_electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(dense_solver, dense_environment).groundStateEnergy();


    // Check if the dense and Davidson energies are equal.
    BOOST_CHECK(std::abs(dense_electronic_energy - davidson_electronic_energy) < 1.0e-08);
}


/**
 *  Check if the ground state energy found using our dense unrestricted FCI routines matches Psi4 and GAMESS' FCI energy.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(unrestricted_FCI_dense) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in a random orthonormal generalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::USpinorBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();

    auto sq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = sq_hamiltonian.numberOfOrbitals();
    sq_hamiltonian.rotate(GQCP::UTransformation<double>::RandomUnitary(K));


    // Set up the full spin-resolved ONV basis.
    GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if the ground state energy found using our dense generalized FCI routines matches Psi4 and GAMESS' FCI energy.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a generalized FCI dimension of 1001, compared to 441 if the correct spin-resolved sector is used.
 */
BOOST_AUTO_TEST_CASE(generalized_FCI_dense) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in a random orthonormal generalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();

    auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(spinor_basis, molecule);
    const auto M = sq_hamiltonian.numberOfOrbitals();
    sq_hamiltonian.rotate(GQCP::GTransformation<double>::RandomUnitary(M));


    // Set up the full spin-unresolved ONV basis.
    GQCP::SpinUnresolvedONVBasis onv_basis {M, molecule.numberOfElectrons()};

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinUnresolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::Operator::NuclearRepulsion(molecule).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Suppose we do a FCI calculation in a random orthonormal basis and we calculate the corresponding 1-DM. Check that, if we rotate to the basis of the FCI naturals, and re-solve the FCI problem, the 1-DM stays the same.
 */
BOOST_AUTO_TEST_CASE(naturals) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.


    // Set up the full spin-resolved ONV basis (with addressing scheme).
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};  // The dimension of this ONV basis is 100.


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto linear_expansion_before = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DM and diagonalize it to obtain the FCI naturals.
    const auto D_before = linear_expansion_before.calculate1DM();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> diagonalizer {D_before};
    GQCP::RTransformation<double> U {diagonalizer.eigenvectors()};


    // Rotate the Hamiltonian to the basis of the FCI naturals, and re-do the FCI calculation. Subsequently check if the 1-DM, calculated from the FCI calculation in the natural orbital basis, is equal to the previously calculated 1-DM's eigenvalues on the diagonal.
    sq_hamiltonian.rotate(U);
    environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    const auto linear_expansion_after = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();

    const auto D_after = linear_expansion_after.calculate1DM();

    BOOST_CHECK(D_after.diagonal().isApprox(diagonalizer.eigenvalues(), 1.0e-12));
}
