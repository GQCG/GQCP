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

#include "WaveFunction/WaveFunctionReader.hpp"


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


BOOST_AUTO_TEST_CASE ( reader_test ) {

    GQCP::WaveFunctionReader test_reader ("../tests/data/test_GAMESS_expansion");
    Eigen::Vector2d test_vector;
    test_vector << 1.0000, 0.0000;

    // Check read vector is correct
    BOOST_CHECK(test_vector.isApprox(test_reader.get_coefficients()));

    // Check if the expansions are equal
    std::string alpha1 = "0000000000000000000000000000000000000000000001";
    std::string alpha2 = "0000000000000000000000000000000000000000000001";
    std::string beta1 = "0000000000000000000000000000000000000000000001";
    std::string beta2 = "0000000000000000000000000000000000000000000010";

    GQCP::Configuration configuration1 = test_reader.get_fock_space().get_Configuration(0);
    GQCP::Configuration configuration2 = test_reader.get_fock_space().get_Configuration(1);

    std::string alpha1_test;
    std::string alpha2_test;
    std::string beta1_test;
    std::string beta2_test;

    std::ostringstream ss;

    ss << configuration1.onv_alpha;
    alpha1_test = ss.str();
    ss.str("");
    ss << configuration1.onv_beta;
    beta1_test = ss.str();
    ss.str("");
    ss << configuration2.onv_alpha;
    alpha2_test = ss.str();
    ss.str("");
    ss << configuration2.onv_beta;
    beta2_test = ss.str();
    ss.str("");

    BOOST_CHECK(alpha1_test == alpha1);
    BOOST_CHECK(alpha2_test == alpha2);
    BOOST_CHECK(beta1_test == beta1);
    BOOST_CHECK(beta2_test == beta2);

}
