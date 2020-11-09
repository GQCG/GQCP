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

#define BOOST_TEST_MODULE "RSpinOrbitalBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check if the orbitals in an AO basis are not orthonormal, but after a transformation to the canonical RHF orbitals, they are.
 */
BOOST_AUTO_TEST_CASE(RHF_orbitals_are_orthonormal) {

    // The orbitals in an AO basis are not orthonormal.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    BOOST_CHECK(!spinor_basis.isOrthonormal());


    // The orbitals in the RHF basis should be orthonormal.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the scalar/AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    spinor_basis.transform(rhf_parameters.coefficientMatrix());
    BOOST_CHECK(spinor_basis.isOrthonormal());
}


/**
 *  Check if a numerical integration over a grid yields overlap values equal to those from Libint. In this test, we won't be using a very high-resolution grid, so so we can't expect very much about the accuracy of the integration.
 * 
 *  The test system is H2//STO-3G.
 */
BOOST_AUTO_TEST_CASE(basisFunctions_integration) {

    // Calculate the overlap matrix through Libint2.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap().parameters();


    // Calculate the basis functions that are in the restricted spin-orbital basis.
    const auto basis_functions = spinor_basis.scalarBasis().basisFunctions();
    const auto& bf1 = basis_functions[0];
    const auto& bf2 = basis_functions[1];


    // Integrate the basis functions over a cubic grid and check if the integration values match Libint's calculations.
    const auto grid = GQCP::CubicGrid::Centered(GQCP::Vector<double, 3>::Zero(), 50, 0.2);

    const auto bf1_squared_evaluated = grid.evaluate(bf1 * bf1);
    const auto bf1_bf2_evaluated = grid.evaluate(bf1 * bf2);
    const auto bf2_squared_evaluated = grid.evaluate(bf2 * bf2);

    BOOST_CHECK(std::abs(grid.integrate(bf1_squared_evaluated) - S(0, 0)) < 1.0e-04);  // 1.0e-04 is still reasonable given the accuracy of the grid
    BOOST_CHECK(std::abs(grid.integrate(bf2_squared_evaluated) - S(1, 1)) < 1.0e-04);
    BOOST_CHECK(std::abs(grid.integrate(bf1_bf2_evaluated) - S(0, 1)) < 1.0e-04);
}


/**
 *  Check if the basis function selector for the Mulliken partitioning generates the correct indices for the AOs.
 */
BOOST_AUTO_TEST_CASE(Mulliken_selector_CO_basis_functions) {

    // Create a test molecule.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const GQCP::Vector<double, 3> C_position {0.0, 0.0, -1.23701};  // In Bohr.
    const GQCP::Vector<double, 3> O_position {0.0, 0.0, 0.927667};  // In Bohr.

    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};


    // Create the reference indices for CO in an STO-3G basis.
    const std::vector<size_t> ref_C_indices {0, 1, 2, 3, 4};  // The indices of the AOs centered on 'C'.
    const std::vector<size_t> ref_O_indices {5, 6, 7, 8, 9};  // The indices of the AOs centered on 'O'.

    // Check if the Mulliken partitioning encapsulates the correct indices.
    using BasisFunction = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>::BasisFunction;
    const auto C_mulliken = spin_orbital_basis.mullikenPartitioning(
        [&C_position](const BasisFunction& bf) {
            return bf.functions()[0].center().isApprox(C_position, 1.0e-04);
        });
    BOOST_REQUIRE(C_mulliken.numberOfOrbitals() == 10);
    BOOST_TEST(C_mulliken.indices() == ref_C_indices, boost::test_tools::per_element());


    const auto O_mulliken = spin_orbital_basis.mullikenPartitioning(
        [&O_position](const BasisFunction& bf) {
            return bf.functions()[0].center().isApprox(O_position, 1.0e-04);
        });
    BOOST_REQUIRE(O_mulliken.numberOfOrbitals() == 10);
    BOOST_TEST(O_mulliken.indices() == ref_O_indices, boost::test_tools::per_element());
}


/**
 *  Check if the shell selector for the Mulliken partitioning generates the correct indices for the AOs.
 */
BOOST_AUTO_TEST_CASE(Mulliken_selector_CO_basis_shells) {

    // Create a test molecule.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};


    // Create the reference indices for CO in an STO-3G basis.
    const std::vector<size_t> ref_C_indices {0, 1, 2, 3, 4};  // The indices of the AOs centered on 'C'.
    const std::vector<size_t> ref_O_indices {5, 6, 7, 8, 9};  // The indices of the AOs centered on 'O'.

    // Check if the Mulliken partitioning encapsulates the correct indices.
    using Shell = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>::Shell;
    const auto C_mulliken = spin_orbital_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "C";
        });
    BOOST_REQUIRE(C_mulliken.numberOfOrbitals() == 10);
    BOOST_TEST(C_mulliken.indices() == ref_C_indices, boost::test_tools::per_element());


    const auto O_mulliken = spin_orbital_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "O";
        });
    BOOST_REQUIRE(O_mulliken.numberOfOrbitals() == 10);
    BOOST_TEST(O_mulliken.indices() == ref_O_indices, boost::test_tools::per_element());
}
