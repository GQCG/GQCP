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
#define BOOST_TEST_MODULE "TransformationMatrix_test"

#include <boost/test/unit_test.hpp>

#include "Basis/TransformationMatrix.hpp"


/**
 *  Check if the 'product' of two transformation matrices behaves as expected
 */
BOOST_AUTO_TEST_CASE ( transform ) {

    // Create two transformation matrices
    GQCP::TransformationMatrix<double> T1 (2);
    T1 << 1.0, 0.0,
          1.0, 3.0;

    GQCP::TransformationMatrix<double> T2 (2);
    T2 << 1.0, 2.0,
          3.0, 4.0;


    // Set up and check the expected result
    GQCP::TransformationMatrix<double> T_total (2);
    T_total <<  1.0,  2.0,
               10.0, 14.0;

    T1.transform(T2);  // first do T1, then T2
    BOOST_CHECK(T1.isApprox(T_total, 1.0e-08));
}
