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
#define BOOST_TEST_MODULE "Hubbard"

#include <boost/test/unit_test.hpp>

#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


/**
 *  Check if a correct constructor does not throw.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_constructor_throws ) {

    const GQCP::SpinResolvedONVBasis onv_basis (15, 3, 3);
    BOOST_CHECK_NO_THROW(GQCP::Hubbard Hubbard {onv_basis});
}


/**
 *  Check if correct arguments for the Hubbard HamiltonianBuilder API throw (and don't) as expected.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_public_methods ) {

    // Create a Hubbard model Hamiltonian.
    const size_t K = 7;  // number of lattice sites
    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


    // Create a compatible ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis_valid (K, 3, 3);
    const auto dim = onv_basis_valid.dimension();
    const GQCP::VectorX<double> x = GQCP::VectorX<double>::Random(dim);

    const GQCP::Hubbard hubbard_builder_valid (onv_basis_valid);
    const auto diagonal = hubbard_builder_valid.calculateDiagonal(hubbard_hamiltonian);


    BOOST_CHECK_NO_THROW(hubbard_builder_valid.constructHamiltonian(hubbard_hamiltonian));
    BOOST_CHECK_NO_THROW(hubbard_builder_valid.matrixVectorProduct(hubbard_hamiltonian, x, diagonal));


    // Create an incompatible ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis_invalid (K+1, 3, 3);

    const GQCP::Hubbard hubbard_builder_invalid (onv_basis_invalid);
    BOOST_CHECK_THROW(hubbard_builder_invalid.constructHamiltonian(hubbard_hamiltonian), std::invalid_argument);
    BOOST_CHECK_THROW(hubbard_builder_invalid.matrixVectorProduct(hubbard_hamiltonian, x, diagonal), std::invalid_argument);
}


/**
 *  Check if the specialized Hubbard Hamiltonian matrix construction algorithm produces the same result as the unspecialized algorithm.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_matrix_specialized_vs_unspecialized ) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 4;  // number of lattice sites
    const size_t N_P = 2;  // number of electron pairs

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian (H);

    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Check if the matrix construction produces identical results.
    const GQCP::Hubbard specialized_builder (onv_basis);
    const auto H_specialized = specialized_builder.constructHamiltonian(hubbard_hamiltonian);

    const GQCP::FCI unspecialized_builder (onv_basis);
    const auto H_unspecialized = unspecialized_builder.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(H_specialized.isApprox(H_unspecialized));
}


/**
 *  Check if the specialized Hubbard Hamiltonian matrix construction algorithm produces the same result as the unspecialized algorithm, for a larger number of lattice sites.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_matrix_specialized_vs_unspecialized_large ) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 6;  // number of lattice sites
    const size_t N_P = 3;  // number of electron pairs

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian (H);

    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Check if the matrix construction produces identical results.
    const GQCP::Hubbard specialized_builder (onv_basis);
    const auto H_specialized = specialized_builder.constructHamiltonian(hubbard_hamiltonian);

    const GQCP::FCI unspecialized_builder (onv_basis);
    const auto H_unspecialized = unspecialized_builder.constructHamiltonian(sq_hamiltonian);

    BOOST_CHECK(H_specialized.isApprox(H_unspecialized));
}


/**
 *  Check if the specialized Hubbard matrix-vector product algorithm produces the same result as the unspecialized algorithm.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_matvec_specialized_vs_unspecialized ) {

        // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 4;  // number of lattice sites
    const size_t N_P = 2;  // number of electron pairs

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian (H);

    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);
    const auto dim = onv_basis.dimension();
    const GQCP::VectorX<double> x = GQCP::VectorX<double>::Random(dim);

    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Check if the matrix-vector product produces identical results.
    const GQCP::Hubbard specialized_builder (onv_basis);
    const auto specialized_diagonal = specialized_builder.calculateDiagonal(hubbard_hamiltonian);
    const auto spezialized_matvec = specialized_builder.matrixVectorProduct(hubbard_hamiltonian, x, specialized_diagonal);

    const GQCP::FCI unspecialized_builder (onv_basis);
    const auto unspecialized_diagonal = unspecialized_builder.calculateDiagonal(sq_hamiltonian);
    const auto unspecialized_matvec = unspecialized_builder.matrixVectorProduct(sq_hamiltonian, x, unspecialized_diagonal);

    BOOST_CHECK(specialized_diagonal.isApprox(unspecialized_diagonal));
    BOOST_CHECK(spezialized_matvec.isApprox(unspecialized_matvec));
}


/**
 *  Check if the specialized Hubbard matrix-vector product algorithm produces the same result as the unspecialized algorithm, for a larger number of lattice sites.
 */
BOOST_AUTO_TEST_CASE ( Hubbard_matvec_specialized_vs_unspecialized_large ) {

        // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 6;  // number of lattice sites
    const size_t N_P = 3;  // number of electron pairs

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian (H);

    GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);
    const auto dim = onv_basis.dimension();
    const GQCP::VectorX<double> x = GQCP::VectorX<double>::Random(dim);

    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Check if the matrix-vector product produces identical results.
    const GQCP::Hubbard specialized_builder (onv_basis);
    const auto specialized_diagonal = specialized_builder.calculateDiagonal(hubbard_hamiltonian);
    const auto spezialized_matvec = specialized_builder.matrixVectorProduct(hubbard_hamiltonian, x, specialized_diagonal);

    const GQCP::FCI unspecialized_builder (onv_basis);
    const auto unspecialized_diagonal = unspecialized_builder.calculateDiagonal(sq_hamiltonian);
    const auto unspecialized_matvec = unspecialized_builder.matrixVectorProduct(sq_hamiltonian, x, unspecialized_diagonal);

    BOOST_CHECK(specialized_diagonal.isApprox(unspecialized_diagonal));
    BOOST_CHECK(spezialized_matvec.isApprox(unspecialized_matvec));
}
