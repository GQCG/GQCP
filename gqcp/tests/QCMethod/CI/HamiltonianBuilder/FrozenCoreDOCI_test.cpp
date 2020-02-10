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
#define BOOST_TEST_MODULE "FrozenCoreDOCI"

#include <boost/test/unit_test.hpp>

#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreDOCI.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/SelectedCI.hpp"


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_constructor ) {

    // Check if a correct constructor works
    GQCP::SpinUnresolvedONVBasis fock_space (15, 3);
    GQCP::SpinUnresolvedFrozenONVBasis frozen_fock_space (fock_space, 2);
    BOOST_CHECK_NO_THROW(GQCP::FrozenCoreDOCI doci (frozen_fock_space));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Random(K);


    // Create a compatible ONV basis
    GQCP::SpinUnresolvedFrozenONVBasis fock_space (K, 3, 1);
    GQCP::FrozenCoreDOCI random_doci (fock_space);
    GQCP::VectorX<double> x = random_doci.calculateDiagonal(sq_hamiltonian);
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian(sq_hamiltonian));
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(sq_hamiltonian, x, x));


    // Create an incompatible ONV basis
    GQCP::SpinUnresolvedFrozenONVBasis fock_space_invalid (K+1, 3, 1);
    GQCP::FrozenCoreDOCI random_doci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_doci_invalid.constructHamiltonian(sq_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(random_doci_invalid.matrixVectorProduct(sq_hamiltonian, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FrozenCoreDOCI ) {

    // Create H-chain HamiltonianParameters to test results
    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (H5, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, H5);  // in an AO basis

    // Create compatible ONV basis
    GQCP::SpinUnresolvedFrozenONVBasis frozen_fock_space (K, 3, 1);
    GQCP::SpinResolvedSelectedONVBasis fock_space (frozen_fock_space);

    // The SpinResolvedSelectedONVBasis includes the same configurations as the SpinUnresolvedFrozenONVBasis
    // These builder instances should return the same results.
    GQCP::SelectedCI random_sci (fock_space);
    GQCP::FrozenCoreDOCI random_frozen_core_doci (frozen_fock_space);

    GQCP::VectorX<double> sci_diagonal = random_sci.calculateDiagonal(sq_hamiltonian);
    GQCP::VectorX<double> doci_diagonal = random_frozen_core_doci.calculateDiagonal(sq_hamiltonian);

    GQCP::VectorX<double> sci_matvec = random_sci.matrixVectorProduct(sq_hamiltonian, sci_diagonal, sci_diagonal);
    GQCP::VectorX<double> doci_matvec = random_frozen_core_doci.matrixVectorProduct(sq_hamiltonian, doci_diagonal, doci_diagonal);

    GQCP::SquareMatrix<double> sci_ham = random_sci.constructHamiltonian(sq_hamiltonian);
    GQCP::SquareMatrix<double> doci_ham = random_frozen_core_doci.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(sci_diagonal.isApprox(doci_diagonal));
    BOOST_CHECK(sci_matvec.isApprox(doci_matvec));
    BOOST_CHECK(sci_ham.isApprox(doci_ham));
}
