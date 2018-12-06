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
#define BOOST_TEST_MODULE "NRDMCalculator_test"


#include "RDM/NRDMCalculator.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( calculateElement ) {

    // Create a test wave function
    size_t M = 3;
    size_t N = 1;
    GQCP::FockSpace fock_space (M, N);

    Eigen::VectorXd coeff (fock_space.get_dimension());
    coeff << 1, 2, -3;


    std::cout << coeff << std::endl;

    // Check some N-RDM values
    GQCP::NRDMCalculator d (fock_space);
    BOOST_CHECK(std::abs(d.calculateElement({0}, {1}, coeff) - 2.0) < 1.0e-12);
}
