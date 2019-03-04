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
#define BOOST_TEST_MODULE "Matrix"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/Matrix.hpp"


BOOST_AUTO_TEST_CASE ( constructor_assignment ) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(3, 3);

    GQCP::MatrixX<double> M1 (A * B);
    GQCP::MatrixX<double> M2 = A + B;
    GQCP::MatrixX<double> M3 = 2*A;
}


BOOST_AUTO_TEST_CASE ( vector_FromFile ) {

    size_t rows = 7;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.dat", rows), std::runtime_error);  // should be 'datA'


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.data", rows));


    // Check that there's an error when the matrix is incompatible with the given file
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("ref_data/h2o_sto-3g_two_electron_horton.data", rows), std::runtime_error);  // can't read in two-electron data in a vector
}



BOOST_AUTO_TEST_CASE ( readVectorFromFile_example ) {

    // Test the read function on a small example
    size_t rows = 4;
    GQCP::VectorX<double> v_ref (rows);
    v_ref << 1.5, -0.2, 0.002, 8.3314;

    auto v = GQCP::VectorX<double>::FromFile("data/small_vector.data", rows);


    BOOST_CHECK(v.isApprox(v_ref, 1.0e-8));
}
