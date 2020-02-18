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
#define BOOST_TEST_MODULE "SelectedCI"

#include <boost/test/unit_test.hpp>

#include "Molecule/Molecule.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/SelectedCI.hpp"


BOOST_AUTO_TEST_CASE ( SelectedCI_constructor ) {

    // Check if a correct constructor works
    GQCP::SpinResolvedONVBasis product_fock_space (5, 3, 3);
    GQCP::SpinResolvedSelectedONVBasis fock_space (product_fock_space);
    BOOST_CHECK_NO_THROW(GQCP::SelectedCI selected_ci (fock_space));
}


BOOST_AUTO_TEST_CASE ( SelectedCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility of the arguments
    size_t K = 5;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Create a compatible ONV basis
    GQCP::SpinResolvedONVBasis product_fock_space (K, 3, 3);
    GQCP::SpinResolvedSelectedONVBasis fock_space (product_fock_space);
    GQCP::SelectedCI random_selected_ci (fock_space);
    GQCP::VectorX<double> x = random_selected_ci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_selected_ci.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_selected_ci.matrixVectorProduct(sq_hamiltonian, x, x));

    // Create an incompatible ONV basis
    GQCP::SpinResolvedONVBasis product_fock_space_invalid (K+1, 3, 3);
    GQCP::SpinResolvedSelectedONVBasis fock_space_invalid (product_fock_space_invalid);
    GQCP::SelectedCI random_selected_ci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_selected_ci_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_selected_ci_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FCI ) {

    // Create H-chain HamiltonianParameters to test results from FCI and selected CI
    size_t K = 4;
    GQCP::Molecule molecule = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in an AO basis

    // Create compatible ONV bases
    GQCP::SpinResolvedONVBasis product_fock_space (K, 2, 2);
    GQCP::SpinResolvedSelectedONVBasis fock_space (product_fock_space);

    // The SpinResolvedSelectedONVBasis includes the same configurations as the SpinResolvedONVBasis, so the HamiltonianBuilders should return the same results
    GQCP::SelectedCI selected_ci (fock_space);
    GQCP::FCI fci (product_fock_space);

    GQCP::VectorX<double> selected_ci_diagonal = selected_ci.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> fci_diagonal = fci.calculateDiagonal(sq_hamiltonian);

    GQCP::VectorX<double> selected_ci_matvec = selected_ci.matrixVectorProduct(sq_hamiltonian, selected_ci_diagonal, selected_ci_diagonal);
    GQCP::VectorX<double> fci_matvec = fci.matrixVectorProduct(sq_hamiltonian, fci_diagonal, fci_diagonal);

    GQCP::SquareMatrix<double> selected_ci_hamiltonian = selected_ci.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> fci_hamiltonian = fci.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(selected_ci_diagonal.isApprox(fci_diagonal));
    BOOST_CHECK(selected_ci_matvec.isApprox(fci_matvec));
    BOOST_CHECK(selected_ci_hamiltonian.isApprox(fci_hamiltonian));
}

/**
 *  Check if DOCI yields the same diagonal, matrix-vector product and Hamiltonian matrix representation as an equivalent selected CI.
 */
BOOST_AUTO_TEST_CASE ( DOCI_vs_selected_CI ) {

    // Set up a Hamiltonian in an orthonormal spinor basis for a linear chain of 4 hydrogens, 1.1 bohr apart.
    const auto molecule = GQCP::Molecule::HChain(4, 1.1);

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Create 'equivalent' ONV bases and 'builders' that 'know' how to generate matrix representations.
    const GQCP::SeniorityZeroONVBasis sz_onv_basis (K, 2);
    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis (sz_onv_basis);

    const GQCP::DOCI doci_builder (sz_onv_basis);
    const GQCP::SelectedCI selected_ci_builder (selected_onv_basis);


    // Check if the calculated diagonals are equal.
    const auto doci_diagonal = doci_builder.calculateDiagonal(sq_hamiltonian);
    const auto selected_ci_diagonal = selected_ci_builder.calculateDiagonal(sq_hamiltonian);

    BOOST_CHECK(doci_diagonal.isApprox(selected_ci_diagonal, 1.0e-12));


    // Check if the matrix-vector products are equal.
    const auto doci_matvec = doci_builder.matrixVectorProduct(sq_hamiltonian, doci_diagonal, doci_diagonal);
    const auto selected_ci_matvec = selected_ci_builder.matrixVectorProduct(sq_hamiltonian, selected_ci_diagonal, selected_ci_diagonal);

    BOOST_CHECK(doci_matvec.isApprox(selected_ci_matvec, 1.0e-12));


    // Check if the matrix representations of the Hamiltonian are equal.
    const auto doci_hamiltonian_matrix = doci_builder.constructHamiltonian(sq_hamiltonian);
    const auto selected_ci_hamiltonian_matrix = selected_ci_builder.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(doci_hamiltonian_matrix.isApprox(selected_ci_hamiltonian_matrix, 1.0e-12));
}
