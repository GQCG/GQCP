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
#define BOOST_TEST_MODULE "CISolver"
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"



BOOST_AUTO_TEST_CASE ( Solver_constructor ) {

    // Create random HamiltonianParameters
    size_t K = 7;
    auto random_hamiltonian_parameters = GQCP::HamiltonianParameters<double>::Random(K);

    // Create a compatible Fock space
    GQCP::FockSpace fock_space (K, 3);
    GQCP::DOCI random_doci (fock_space);
    BOOST_CHECK_NO_THROW(GQCP::CISolver ci_solver (random_doci, random_hamiltonian_parameters));

    // Create an incompatible Fock space
    GQCP::FockSpace fock_space_invalid (K+1, 3);
    GQCP::DOCI random_doci_invalid (fock_space_invalid);
    BOOST_CHECK_THROW(GQCP::CISolver ci_solver (random_doci_invalid, random_hamiltonian_parameters), std::invalid_argument);
}
