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
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "HamiltonianBuilder/FrozenCoreDOCI.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/SelectedCI.hpp"




BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_constructor ) {

    // Check if a correct constructor works
    GQCP::FockSpace fock_space (15, 3);
    GQCP::FrozenFockSpace frozen_fock_space (fock_space, 2);
    BOOST_CHECK_NO_THROW(GQCP::FrozenCoreDOCI doci (frozen_fock_space));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_public_methods ) {

    // Create random HamiltonianParameters to check compatibility
    size_t K = 5;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Random(K);


    // Create a compatible Fock space
    GQCP::FrozenFockSpace fock_space (K, 3, 1);
    GQCP::FrozenCoreDOCI random_doci (fock_space);
    Eigen::VectorXd x = random_doci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(random_hamiltonian_parameters, x, x));


    // Create an incompatible Fock space
    GQCP::FrozenFockSpace fock_space_invalid (K+1, 3, 1);
    GQCP::FrozenCoreDOCI random_doci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(random_doci_invalid.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_doci_invalid.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( SelectedCI_vs_FrozenCoreDOCI ) {

    // Create H-chain HamiltonianParameters to test results
    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    // Create compatible Fock spaces
    GQCP::FrozenFockSpace frozen_fock_space (K, 3, 1);
    GQCP::SelectedFockSpace fock_space (frozen_fock_space);

    // The SelectedFockSpace includes the same configurations as the FrozenFockSpace
    // These builder instances should return the same results.
    GQCP::SelectedCI random_sci (fock_space);
    GQCP::FrozenCoreDOCI random_frozen_core_doci (frozen_fock_space);

    Eigen::VectorXd sci_diagonal = random_sci.calculateDiagonal(random_hamiltonian_parameters);
    Eigen::VectorXd doci_diagonal = random_frozen_core_doci.calculateDiagonal(random_hamiltonian_parameters);

    Eigen::VectorXd sci_matvec = random_sci.matrixVectorProduct(random_hamiltonian_parameters, sci_diagonal, sci_diagonal);
    Eigen::VectorXd doci_matvec = random_frozen_core_doci.matrixVectorProduct(random_hamiltonian_parameters, doci_diagonal, doci_diagonal);

    Eigen::MatrixXd sci_ham = random_sci.constructHamiltonian(random_hamiltonian_parameters);
    Eigen::MatrixXd doci_ham = random_frozen_core_doci.constructHamiltonian(random_hamiltonian_parameters);

    BOOST_CHECK(sci_diagonal.isApprox(doci_diagonal));
    BOOST_CHECK(sci_matvec.isApprox(doci_matvec));
    BOOST_CHECK(sci_ham.isApprox(doci_ham));
}
