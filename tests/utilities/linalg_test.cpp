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
#define BOOST_TEST_MODULE "linalg_test"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "utilities/linalg.hpp"


BOOST_AUTO_TEST_CASE ( areEqualEigenvectors ) {

    // Test areEqualEigenvectors with an example.
    GQCP::VectorX<double> a (3);
    GQCP::VectorX<double> b (3);
    GQCP::VectorX<double> c (3);
    GQCP::VectorX<double> d (3);

    a << 2, 3, 1;
    b << 2, 3, 1;
    c << -2, -3, -1;
    d << 2, 3, 0;


    BOOST_CHECK(GQCP::areEqualEigenvectors(a, b, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualEigenvectors(a, c, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualEigenvectors(b, c, 1.0e-6));

    BOOST_CHECK(!GQCP::areEqualEigenvectors(a, d, 1.0e-6));
    BOOST_CHECK(!GQCP::areEqualEigenvectors(c, d, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqualSetsOfEigenvectors_throws ) {

    // Check for throws if the dimensions aren't compatible.
    GQCP::MatrixX<double> C1 (3, 3);
    GQCP::MatrixX<double> C2 (3, 2);

    BOOST_CHECK_THROW(GQCP::areEqualSetsOfEigenvectors(C1, C2, 1.0e-6), std::invalid_argument);


    // Check for no throw if the dimensions are compatible
    GQCP::MatrixX<double> C3 (3, 3);

    BOOST_CHECK_NO_THROW(GQCP::areEqualSetsOfEigenvectors(C1, C3, 1.0e-6));
}


BOOST_AUTO_TEST_CASE ( areEqualSetsOfEigenvectors_example ) {

    // Test areEqualSetsOfEigenvectors with an example
    GQCP::MatrixX<double> eigenvectors1 (2, 2);
    GQCP::MatrixX<double> eigenvectors2 (2, 2);
    GQCP::MatrixX<double> eigenvectors3 (2, 2);
    GQCP::MatrixX<double> eigenvectors4 (2, 2);

    eigenvectors1 << 0,  2,
                     1, -1;

    eigenvectors2 << 0,  2,
                     1, -1;

    eigenvectors3 << 0, -2,
                     1,  1;

    eigenvectors4 << 1,  2,
                     1, -1;


    BOOST_CHECK(GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors2, 1.0e-6));
    BOOST_CHECK(GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors3, 1.0e-6));

    BOOST_CHECK(!GQCP::areEqualSetsOfEigenvectors(eigenvectors1, eigenvectors4, 1.0e-6));
}
