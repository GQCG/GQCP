// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "Field_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Grid/Field.hpp"


/**
 *  Check if reading in the data in a GAUSSIAN Cube file is correct.
 */
BOOST_AUTO_TEST_CASE(ReadCubeFile) {

    // Read in the Cube file, provide the reference values and check if it was parsed correctly.
    const auto scalar_field = GQCP::Field<double>::ReadCubeFile("data/benzene.cube");


    // Check the results.
    BOOST_CHECK(scalar_field.size() == 60 * 60 * 60);
    BOOST_CHECK(std::abs(scalar_field.value(0) - (-1.38100446455e-06)) < 1.0e-08);
    BOOST_CHECK(std::abs(scalar_field.values().back() - (-9.58547463676e-07)) < 1.0e-08);
}


/**
 *  Check if reading in an integration grid file works as expected.
 */
BOOST_AUTO_TEST_CASE(ReadIntegrationGridFile) {

    // Read in the integration grid file and check if the parsing went correctly.
    const auto J = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/benzene_current.igrid");

    BOOST_CHECK(J.size() == 25905);


    // Check the first field value.
    GQCP::Vector<double, 3> ref_value_first {0, 0, 0};
    BOOST_CHECK(J.value(0).isApprox(ref_value_first, 1.0e-08));

    // Check the last field value.
    GQCP::Vector<double, 3> ref_value_last {-1.4694209851333318E-003, 5.1869592466737078E-004, 1.2283072801926104E-005};
    BOOST_CHECK(J.values().back().isApprox(ref_value_last, 1.0e-08));
}


/**
 *  Check if reading in a regular grid file works as expected.
 */
BOOST_AUTO_TEST_CASE(ReadRegularGridFile) {

    // Read in the regular grid file and check if the parsing went correctly.
    const auto J = GQCP::Field<GQCP::Vector<double, 3>>::ReadGridFile<3>("data/benzene_current.rgrid");

    BOOST_CHECK(J.size() == 13824);


    // Check the first field value.
    GQCP::Vector<double, 3> ref_value_first {-0.000002, 0.000002, 0.000000};
    BOOST_CHECK(J.value(0).isApprox(ref_value_first, 1.0e-08));

    // Check the last field value.
    GQCP::Vector<double, 3> ref_value_last {0.000003, -0.000004, -0.000001};
    BOOST_CHECK(J.values().back().isApprox(ref_value_last, 1.0e-08));
}


/**
 *  Check if operator+ is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(operator_plus) {

    // Set up two fields.
    GQCP::Field<unsigned> lhs {{8, 5, 5}};
    GQCP::Field<unsigned> rhs {{1, 2, 3}};


    // Check if operator+ works as expected.
    const std::vector<unsigned> ref_sum_values1 {9, 7, 8};
    auto sum = lhs + rhs;
    BOOST_CHECK_EQUAL_COLLECTIONS(sum.values().begin(), sum.values().end(), ref_sum_values1.begin(), ref_sum_values1.end());


    // Check if operator+= works as expected.
    const std::vector<unsigned> ref_sum_values2 {10, 9, 11};
    sum += rhs;
    BOOST_CHECK_EQUAL_COLLECTIONS(sum.values().begin(), sum.values().end(), ref_sum_values2.begin(), ref_sum_values2.end());
}


/**
 *  Check if operator- is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(operator_minus) {

    // Set up two fields.
    GQCP::Field<int> lhs {{8, 5, 5}};
    GQCP::Field<int> rhs {{1, 2, 3}};


    // Check if the unary operator- works as expected.
    const std::vector<int> ref_unary_minus_values {-8, -5, -5};
    const auto unary_minus = -lhs;
    BOOST_CHECK_EQUAL_COLLECTIONS(unary_minus.values().begin(), unary_minus.values().end(), ref_unary_minus_values.begin(), ref_unary_minus_values.end());


    // Check if binary operator- works as expected.
    const std::vector<int> ref_difference_values1 {7, 3, 2};
    auto difference = lhs - rhs;
    BOOST_CHECK_EQUAL_COLLECTIONS(difference.values().begin(), difference.values().end(), ref_difference_values1.begin(), ref_difference_values1.end());


    // Check if operator-= works as expected.
    const std::vector<int> ref_difference_values2 {6, 1, -1};
    difference -= rhs;
    BOOST_CHECK_EQUAL_COLLECTIONS(difference.values().begin(), difference.values().end(), ref_difference_values2.begin(), ref_difference_values2.end());
}
