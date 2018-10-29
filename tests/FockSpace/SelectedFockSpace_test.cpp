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

    GQCP::ProductFockSpace fock_space_product (10, 5, 5);
    GQCP::FockSpace fock_space (10, 5);

    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace fock (fock_space_product));
    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace fock (fock_space));
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

    fock_space.addConfiguration(alpha_set, beta_set);


    // Check if the expansions are equal
    // Generate the expected results
    std::string alpha1_ref = "001";
    std::string alpha2_ref = "010";
    std::string beta1_ref = "001";
    std::string beta2_ref = "010";

    // Retrieve the added results
    GQCP::Configuration configuration1 = fock_space.get_configuration(0);
    GQCP::Configuration configuration2 = fock_space.get_configuration(1);

    // Retrieve the string representation of the ONVs
    std::string alpha1_test = configuration1.onv_alpha.asString();
    std::string alpha2_test = configuration2.onv_alpha.asString();
    std::string beta1_test = configuration1.onv_beta.asString();
    std::string beta2_test = configuration2.onv_beta.asString();

    BOOST_CHECK(alpha1_test == alpha1_ref);
    BOOST_CHECK(alpha2_test == alpha2_ref);
    BOOST_CHECK(beta1_test == beta1_ref);
    BOOST_CHECK(beta2_test == beta2_ref);
}


BOOST_AUTO_TEST_CASE ( reader_test ) {

    // We will test if we can construct a selected fock space and a corresponding coefficients
    // through an input file
    GQCP::WaveFunctionReader test_reader ("../tests/data/test_GAMESS_expansion");


    // Check read vector is correct
    // Gerenate the expected result
    Eigen::Vector2d test_vector;
    test_vector << 1.0000, 0.0000;

    BOOST_CHECK(test_vector.isApprox(test_reader.get_coefficients()));

    // Check if the expansions are equal
    // Generate the expected results
    std::string alpha1_ref = "0000000000000000000000000000000000000000000001";
    std::string alpha2_ref = "0000000000000000000000000000000000000000000001";
    std::string beta1_ref = "0000000000000000000000000000000000000000000001";
    std::string beta2_ref = "0000000000000000000000000000000000000000000010";

    // Retrieve the read results
    GQCP::Configuration configuration1 = test_reader.get_fock_space().get_configuration(0);
    GQCP::Configuration configuration2 = test_reader.get_fock_space().get_configuration(1);

    // Retrieve the string representation of the ONVs
    std::string alpha1_test = configuration1.onv_alpha.asString();
    std::string alpha2_test = configuration2.onv_alpha.asString();
    std::string beta1_test = configuration1.onv_beta.asString();
    std::string beta2_test = configuration2.onv_beta.asString();

    BOOST_CHECK(alpha1_test == alpha1_ref);
    BOOST_CHECK(alpha2_test == alpha2_ref);
    BOOST_CHECK(beta1_test == beta1_ref);
    BOOST_CHECK(beta2_test == beta2_ref);

}
