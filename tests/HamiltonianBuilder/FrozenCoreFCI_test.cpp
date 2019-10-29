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
#define BOOST_TEST_MODULE "FrozenCoreFCI"

#include <boost/test/unit_test.hpp>

#include "HamiltonianBuilder/FrozenCoreFCI.hpp"
#include "HamiltonianBuilder/SelectedCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


BOOST_AUTO_TEST_CASE ( FrozenCoreFCI_constructor ) {

    // Check if a correct constructor works
    GQCP::ProductFockSpace fock_space (15, 3, 3);
    GQCP::FrozenProductFockSpace frozen_fock_space (fock_space, 2);
    BOOST_CHECK_NO_THROW(GQCP::FrozenCoreFCI fci (frozen_fock_space));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreFCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);


    // Create a compatible Fock space
    GQCP::FrozenProductFockSpace fock_space (K, 3, 3, 1);
    GQCP::FrozenCoreFCI random_fci (fock_space);
    GQCP::VectorX<double> x = random_fci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_fci.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_fci.matrixVectorProduct(sq_hamiltonian, x, x));


    // Create an incompatible Fock space
    GQCP::FrozenProductFockSpace fock_space_invalid (K+1, 3, 3, 1);
    GQCP::FrozenCoreFCI random_fci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_fci_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_fci_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FrozenCoreFCI ) {

    // Create H-chain HamiltonianParameters to test results
    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (H5, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H5);  // in an AO basis

    // Create compatible Fock spaces
    GQCP::FrozenProductFockSpace product_fock_space (K, 3, 3, 1);
    GQCP::SelectedFockSpace fock_space (product_fock_space);

    // The SelectedFockSpace includes the same configurations as the FrozenProductFockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI random_sci (fock_space);
    GQCP::FrozenCoreFCI random_frozen_core_fci (product_fock_space);

    GQCP::VectorX<double> sci_diagonal = random_sci.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> fci_diagonal = random_frozen_core_fci.calculateDiagonal(sq_hamiltonian);

    GQCP::VectorX<double> sci_matvec = random_sci.matrixVectorProduct(sq_hamiltonian, sci_diagonal, sci_diagonal);
    GQCP::VectorX<double> fci_matvec = random_frozen_core_fci.matrixVectorProduct(sq_hamiltonian, fci_diagonal, fci_diagonal);

    GQCP::SquareMatrix<double> sci_ham = random_sci.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> fci_ham = random_frozen_core_fci.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(sci_diagonal.isApprox(fci_diagonal));
    BOOST_CHECK(sci_matvec.isApprox(fci_matvec));
    BOOST_CHECK(sci_ham.isApprox(fci_ham));
}
