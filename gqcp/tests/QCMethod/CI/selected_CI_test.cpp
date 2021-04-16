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

#define BOOST_TEST_MODULE "SpinResolvedSelectedONVBasis"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCMethod/HF/UHF/UHF.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

/**
 *  Check if the ground state energy found using our restricted selected FCI routines matches Psi4 and GAMESS' FCI energy.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(restricted_selected_FCI) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto K = sq_hamiltonian.numberOfOrbitals();

    // Set up the full spin-resolved selected ONV basis.
    GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, selected_onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}


/**
 *  Check if the ground state energy found using our unrestricted selected FCI routines matches Psi4 and GAMESS' FCI energy.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(unrestricted_selected_FCI) {

    const double reference_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in a random orthonormal unrestriced spin-orbital basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();

    auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto K = sq_hamiltonian.numberOfOrbitals();
    sq_hamiltonian.rotate(GQCP::UTransformation<double>::RandomUnitary(K));


    // Set up the full spin-resolved selected ONV basis.
    GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, selected_onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedSelectedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}

#include <iomanip>

/**
 */
BOOST_AUTO_TEST_CASE(H2O_CIS) {

    std::cout << std::setprecision(15) << std::endl;

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_cis.xyz");

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    // spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    // const auto K = hamiltonian.numberOfOrbitals();

    // Solve the RHF SCF equations to find an initial orthonormal basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    std::cout << "SCF ENERGY: " << rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value() << std::endl;
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    hamiltonian.transform(rhf_parameters.expansion());


    // Set up the full spin-resolved selected ONV basis.
    const auto onv_basis = GQCP::SpinResolvedSelectedONVBasis::CIS(7, 5, 5);


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto ci_qc_structure = GQCP::QCMethod::CI<GQCP::SpinResolvedSelectedONVBasis>(onv_basis, onv_basis.dimension()).optimize(solver, environment);

    for (size_t i = 0; i < onv_basis.dimension(); i++) {
        std::cout << ci_qc_structure.energy(i) + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value() << std::endl;
    }


    // Check our result with the reference.
    // BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}
