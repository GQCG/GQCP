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
#define BOOST_TEST_MODULE "ChemicalMatrix"

#include <boost/test/unit_test.hpp>

#include "Mathematical/ChemicalMatrix.hpp"

#include <boost/math/constants/constants.hpp>


BOOST_AUTO_TEST_CASE ( constructor ) {

    auto M1 = GQCP::ChemicalMatrix<double>::Zero(2, 2);
    BOOST_CHECK_NO_THROW(GQCP::ChemicalMatrix<double> square_M1 (M1));

    auto M2 = GQCP::ChemicalMatrix<double>::Zero(2, 1);
    BOOST_CHECK_THROW(GQCP::ChemicalMatrix<double> square_M2 (M2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( ChemicalMatrix_transform_trivial ) {

    // Let's test a trivial transformation: i.e. with T being a unit matrix
    GQCP::ChemicalMatrix<double> h = GQCP::ChemicalMatrix<double>::Random(3, 3);
    GQCP::ChemicalMatrix<double> H = h;

    GQCP::ChemicalMatrix<double> T = GQCP::ChemicalMatrix<double>::Identity(3, 3);
    H.basisTransform(T);

    BOOST_CHECK(H.isApprox(h, 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( ChemicalMatrix_transform_and_inverse ) {

    // Let's test if, if we transform h with T and then with T_inverse, we get effectively do nothing
    GQCP::ChemicalMatrix<double> h = GQCP::ChemicalMatrix<double>::Random(3, 3);
    GQCP::ChemicalMatrix<double> H = h;

    GQCP::ChemicalMatrix<double> T (3);
    T <<    1,  0,  0,
            0, -2,  0,
            0,  0,  3;
    GQCP::ChemicalMatrix<double> T_inverse = T.inverse();


    H.basisTransform(T);
    H.basisTransform(T_inverse);

    BOOST_CHECK(H.isApprox(h, 1.0e-12));
}
