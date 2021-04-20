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
    auto solver = GQCP::EigenproblemSolver::Dense<double>();
    const auto electronic_energy = GQCP::QCMethod::CI<double, GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

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
    auto solver = GQCP::EigenproblemSolver::Dense<double>();
    const auto electronic_energy = GQCP::QCMethod::CI<double, GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();

    // Check our result with the reference.
    const auto energy = electronic_energy + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-06);
}
