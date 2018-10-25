// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "FCI"


#include "HamiltonianBuilder/FCI.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( FCI_constructor ) {
    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (15, 3, 3);

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCP::FCI fci (fock_space));
}


BOOST_AUTO_TEST_CASE ( FCI_public_methods ) {
    // Create an AOBasis
    GQCP::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis->get_number_of_basis_functions();
    GQCP::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCP::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCP::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCP::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);

    // Create a compatible Fock space
    GQCP::ProductFockSpace fock_space (K, 3, 3);

    // Create FCI module
    GQCP::FCI random_fci (fock_space);

    // Test the public FCI methods
    Eigen::VectorXd x = random_fci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_fci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_fci.matrixVectorProduct(random_hamiltonian_parameters, x, x));

    // Create an incompatible Fock space
    GQCP::ProductFockSpace fock_space_i (K+1, 3, 3);

    // Create FCI module
    GQCP::FCI random_fci_i (fock_space_i);
    BOOST_CHECK_THROW(random_fci_i.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_fci_i.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}
