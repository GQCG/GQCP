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


/**
 *  Check if TransformationMatrix can correctly express itself in the GAMESS-US $VEC group notation 
 */
BOOST_AUTO_TEST_CASE ( gamessUS ) {

	std::string reference_T1 = "$VEC\n1  1 1.00000000E+00 1.00000000E+00\n2  1 0.00000000E+00 3.00000000E+00\n$END\n";
	std::string reference_T2 = "$VEC\n1  1-1.00000000E+00 1.00000000E+00 0.00000000E+00 1.00000000E+00 0.00000000E+00\n1  2 1.00000000E+00\n2  1 0.00000000E+00 3.00000000E+00 0.00000000E+00 9.00000000E+00 0.00000000E+00\n2  2-7.00000000E+00\n3  1-4.00000000E+00-4.00000000E+00-4.00000000E+00-4.00000000E+00-4.00000000E+00\n3  2-4.00000000E+00\n4  1 5.00000000E+00 5.00000000E+00 5.00000000E+00 5.00000000E+00 5.00000000E+00\n4  2 5.00000000E+00\n5  1-7.00000000E+00-7.00000000E+00-7.00000000E+00-7.00000000E+00-7.00000000E+00\n5  2-7.00000000E+00\n6  1 8.00000000E+00 8.00000000E+00 8.00000000E+00 8.00000000E+00 8.00000000E+00\n6  2 8.00000000E+00\n$END\n";

    // Create transformation matrices
    GQCP::TransformationMatrix<double> T1 (2);
    T1 << 1.0, 0.0,
          1.0, 3.0;

    GQCP::TransformationMatrix<double> T2 (6);
	T2 <<-1.0, 0.0,-4.0, 5.0,-7.0, 8.0,
		  1.0, 3.0,-4.0, 5.0,-7.0, 8.0,
		  0.0, 0.0,-4.0, 5.0,-7.0, 8.0,
		  1.0, 9.0,-4.0, 5.0,-7.0, 8.0,
		  0.0, 0.0,-4.0, 5.0,-7.0, 8.0,
		  1.0,-7.0,-4.0, 5.0,-7.0, 8.0;
		
	// Test them against a reference representation
    BOOST_CHECK(reference_T1.compare(T1.asGAMESSUSVECGroup()) == 0);
    BOOST_CHECK(reference_T2.compare(T2.asGAMESSUSVECGroup()) == 0);
}
