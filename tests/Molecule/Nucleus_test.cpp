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
#define BOOST_TEST_MODULE "Nucleus"

#include <boost/test/unit_test.hpp>

#include "Molecule/Nucleus.hpp"



BOOST_AUTO_TEST_CASE ( Nucleus_constructor ) {

    GQCP::Nucleus nucleus {1, 0.0, 0.1, 0.2};
}


BOOST_AUTO_TEST_CASE ( Nucleus_isSmallerThan ) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {2, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus3 {2, 0.1, 0.2, 0.2};
    GQCP::Nucleus nucleus4 {2, 0.1, 0.2, 0.3};


    // Check if operator< does what is expected
    BOOST_CHECK(nucleus1.isSmallerThan(nucleus2));
    BOOST_CHECK(nucleus2.isSmallerThan(nucleus3));
    BOOST_CHECK(nucleus3.isSmallerThan(nucleus4));

    BOOST_CHECK(!(nucleus2.isSmallerThan(nucleus1)));
    BOOST_CHECK(!(nucleus3.isSmallerThan(nucleus2)));
    BOOST_CHECK(!(nucleus4.isSmallerThan(nucleus3)));


    // Check if the tolerance works
    BOOST_CHECK(!(nucleus2.isSmallerThan(nucleus3, 0.2)));
    BOOST_CHECK(!(nucleus3.isSmallerThan(nucleus2, 0.2)));
}


BOOST_AUTO_TEST_CASE ( Nucleus_operator_smaller_than ) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {2, 0.0, 0.1, 0.2};

    // A small test to check if we can operator<
    BOOST_CHECK(nucleus1 < nucleus2);
}


BOOST_AUTO_TEST_CASE ( Nucleus_isEqualTo ) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus3 {2, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus4 {1, 0.1, 0.2, 0.3};

    // Check if they're equal
    BOOST_CHECK(nucleus1.isEqualTo(nucleus2));

    // Check if different nucleusic numbers cause inequality
    BOOST_CHECK(!(nucleus1.isEqualTo(nucleus3)));

    // Check if different coordinates cause inequality
    BOOST_CHECK(!(nucleus1.isEqualTo(nucleus4)));


    // Check if the tolerance works
    BOOST_CHECK(nucleus1.isEqualTo(nucleus2, 0.2));
}


BOOST_AUTO_TEST_CASE ( Nucleus_operator_equals ) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {1, 0.0, 0.1, 0.2};

    // A small test to check if we can operator==
    BOOST_CHECK(nucleus1 == nucleus2);
}


BOOST_AUTO_TEST_CASE ( Nucleus_operator_ostream ) {

    GQCP::Nucleus nucleus1 {1, 0.0, 0.1, 0.2};
    GQCP::Nucleus nucleus2 {2, 0.1, 0.2, 0.3};


    std::cout << nucleus1 << std::endl;
    std::cout << nucleus2 << std::endl;
}


BOOST_AUTO_TEST_CASE ( calculateDistance ) {

    // Create some nuclei
    GQCP::Nucleus nucleus1 {1, 0, 3, 0};
    GQCP::Nucleus nucleus2 {1, 0, 0, 4};
    GQCP::Nucleus nucleus3 {1, 3, 0, 0};
    GQCP::Nucleus nucleus4 {1, 0, 0, 5};


    // Check their distances
    BOOST_CHECK(std::abs(nucleus1.calculateDistance(nucleus2) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistance(nucleus3) - std::sqrt(18.0)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistance(nucleus2) - nucleus2.calculateDistance(nucleus3)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistance(nucleus3) - 5) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistance(nucleus4) - 1) < 1.0e-12);

    // Check that the distances are symmetric
    BOOST_CHECK(std::abs(nucleus1.calculateDistance(nucleus2) - nucleus2.calculateDistance(nucleus1)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus1.calculateDistance(nucleus3) - nucleus3.calculateDistance(nucleus1)) < 1.0e-12);
    BOOST_CHECK(std::abs(nucleus2.calculateDistance(nucleus3) - nucleus3.calculateDistance(nucleus2)) < 1.0e-12);
}
