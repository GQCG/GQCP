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
#define BOOST_TEST_MODULE "SelectedFockSpace"


#include "FockSpace/SelectedFockSpace.hpp"

//#include "WaveFunction/WaveFunctionReader.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace (10, 5, 5));
}


BOOST_AUTO_TEST_CASE ( addConfiguration ) {

    // Create a faulty expansion: one of the orbitals is different
    GQCP::SelectedFockSpace fock_space (3, 1, 1);

    std::vector<std::string> alpha_set = {"001", "010"};
    std::vector<std::string> beta_set = {"001", "010"};

    BOOST_CHECK_NO_THROW(fock_space.addConfiguration(alpha_set, beta_set));

    // Test throw with one of the sets is not the same size
    std::vector<std::string> beta_set_long = {"001", "010", "100"};
    BOOST_CHECK_THROW(fock_space.addConfiguration(alpha_set, beta_set_long), std::invalid_argument);

    // Test throw with incompatible orbital numbers
    BOOST_CHECK_THROW(fock_space.addConfiguration("0001", "0100"), std::invalid_argument);

    // Test throw with incompatible electron numbers
    BOOST_CHECK_THROW(fock_space.addConfiguration("011", "011"), std::invalid_argument);
}

/*
BOOST_AUTO_TEST_CASE ( reader_test ) {

    GQCP::WaveFunctionReader test_reader ("../tests/data/test_GAMESS_expansion");
    Eigen::VectorXd test = {1.0000, 0.0000};

    // Check if both expansions are considered equal
    BOOST_CHECK(test.isApprox(test_reader.get_coefficients()));
}
 */