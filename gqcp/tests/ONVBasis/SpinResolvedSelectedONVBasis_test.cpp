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

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


/**
 *  Test the general functionality of the `expandWith` function, by testing throws and retrieving ONVs.
 */
BOOST_AUTO_TEST_CASE(expandWith) {

    // Set up an initial zero ONV basis.
    GQCP::SpinResolvedSelectedONVBasis onv_basis {3, 1, 1};


    // Add two compatible ONVs.
    BOOST_CHECK_NO_THROW(onv_basis.expandWith(GQCP::SpinResolvedONV::FromString("001", "001")));
    BOOST_CHECK_NO_THROW(onv_basis.expandWith(GQCP::SpinResolvedONV::FromString("010", "010")));


    // Check if we receive a throw if we add an ONV with an incompatible number of orbitals.
    BOOST_CHECK_THROW(onv_basis.expandWith(GQCP::SpinResolvedONV::FromString("0001", "0100")), std::invalid_argument);


    // Check if we receive a throw if we add an ONV with an incompatible number of electrons.
    BOOST_CHECK_THROW(onv_basis.expandWith(GQCP::SpinResolvedONV::FromString("011", "011")), std::invalid_argument);


    // Check if we can access the ONVs in order.
    BOOST_CHECK(onv_basis.onvWithIndex(0).asString() == "001|001");
    BOOST_CHECK(onv_basis.onvWithIndex(1).asString() == "010|010");
}


/**
 *  Check if the matrix-vector product through a direct evaluation (i.e. through the dense Hamiltonian matrix representation) and the specialized implementation are equal.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(restricted_dense_vs_matvec) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-resolved selected ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};

    // Determine the Hamiltonian matrix and let it act on a random linear expansion.
    const auto linear_expansion = GQCP::LinearExpansion<double, GQCP::SpinResolvedSelectedONVBasis>::Random(selected_onv_basis);
    const auto H_dense = selected_onv_basis.evaluateOperatorDense(hamiltonian);
    const GQCP::VectorX<double> direct_mvp = H_dense * linear_expansion.coefficients();  // mvp: matrix-vector-product

    // Determine the specialized matrix-vector product and check if they are equal.
    const auto specialized_mvp = selected_onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, linear_expansion.coefficients());

    BOOST_CHECK(specialized_mvp.isApprox(direct_mvp, 1.0e-08));
}


/**
 *  Check if the diagonal of the matrix representation of a restricted Hamiltonian is equal to the diagonal that is calculated through a specialized routine.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(restricted_hamiltonian_diagonal) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-resolved selected ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};

    // Determine the Hamiltonian matrix and the diagonal through a specialized routine, and check if they match.
    const auto dense_matrix = selected_onv_basis.evaluateOperatorDense(hamiltonian);
    const auto diagonal_specialized = selected_onv_basis.evaluateOperatorDiagonal(hamiltonian);

    BOOST_CHECK(diagonal_specialized.isApprox(dense_matrix.diagonal(), 1.0e-08));
}
