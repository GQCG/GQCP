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

#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/SelectedCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Molecule/Molecule.hpp"


BOOST_AUTO_TEST_CASE ( SelectedCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace product_fock_space (5, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space);
    BOOST_CHECK_NO_THROW(GQCP::SelectedCI selected_ci (fock_space));
}


BOOST_AUTO_TEST_CASE ( SelectedCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility of the arguments
    size_t K = 5;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);

    // Create a compatible Fock space
    GQCP::ProductFockSpace product_fock_space (K, 3, 3);
    GQCP::SelectedFockSpace fock_space (product_fock_space);
    GQCP::SelectedCI random_selected_ci (fock_space);
    GQCP::VectorX<double> x = random_selected_ci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_selected_ci.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_selected_ci.matrixVectorProduct(sq_hamiltonian, x, x));

    // Create an incompatible Fock space
    GQCP::ProductFockSpace product_fock_space_invalid (K+1, 3, 3);
    GQCP::SelectedFockSpace fock_space_invalid (product_fock_space_invalid);
    GQCP::SelectedCI random_selected_ci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_selected_ci_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_selected_ci_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FCI ) {

    // Create H-chain HamiltonianParameters to test results from FCI and selected CI
    size_t K = 4;
    GQCP::Molecule H4 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (H4, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H4);  // in an AO basis

    // Create compatible Fock spaces
    GQCP::ProductFockSpace product_fock_space (K, 2, 2);
    GQCP::SelectedFockSpace fock_space (product_fock_space);

    // The SelectedFockSpace includes the same configurations as the ProductFockSpace, so the HamiltonianBuilders should return the same results
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

BOOST_AUTO_TEST_CASE ( SelectedCI_vs_DOCI ) {

    // Create H-chain HamiltonianParameters to test results from DOCI and selected CI
    size_t K = 4;
    GQCP::Molecule H4 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (H4, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H4);  // in an AO basis

    // Create compatible Fock spaces
    GQCP::FockSpace do_fock_space (K, 2);
    GQCP::SelectedFockSpace fock_space (do_fock_space);

    // The SelectedFockSpace includes the same configurations as the FockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI selected_ci (fock_space);
    GQCP::DOCI doci (do_fock_space);

    GQCP::VectorX<double> selected_ci_diagonal = selected_ci.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> doci_diagonal = doci.calculateDiagonal(sq_hamiltonian);

    GQCP::VectorX<double> selected_ci_matvec = selected_ci.matrixVectorProduct(sq_hamiltonian, selected_ci_diagonal, selected_ci_diagonal);
    GQCP::VectorX<double> doci_matvec = doci.matrixVectorProduct(sq_hamiltonian, doci_diagonal, doci_diagonal);

    GQCP::SquareMatrix<double> selected_ci_hamiltonian = selected_ci.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> doci_hamiltonian = doci.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(selected_ci_diagonal.isApprox(doci_diagonal));
    BOOST_CHECK(selected_ci_matvec.isApprox(doci_matvec));
    BOOST_CHECK(selected_ci_hamiltonian.isApprox(doci_hamiltonian));
}
