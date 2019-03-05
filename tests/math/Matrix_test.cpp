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


BOOST_AUTO_TEST_CASE ( Vector_FromFile ) {

    size_t rows = 4;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.dat", rows), std::runtime_error);  // should be 'datA'


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.data", rows));


    // Check that there's an error when trying to read in tensor data into a vector
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("ref_data/h2o_sto-3g_two_electron_horton.data", rows), std::runtime_error);


    // Test the read function on a small example
    GQCP::VectorX<double> v_ref (rows);
    v_ref << 1.5, -0.2, 0.002, 8.3314;

    auto v = GQCP::VectorX<double>::FromFile("data/small_vector.data", rows);

    BOOST_CHECK(v.isApprox(v_ref, 1.0e-15));
}


BOOST_AUTO_TEST_CASE ( Matrix_FromFile ) {

    size_t rows = 2;
    size_t cols = 2;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::Matrix<double>::FromFile("data/h2o_sto-3g_kinetic_horton.dat", rows, cols), std::runtime_error);  // should be 'datA'


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::Matrix<double>::FromFile("data/small_one_ints.data", rows, cols));


    // Check that there's an error when trying to read in tensor data into a matrix
    BOOST_CHECK_THROW(GQCP::Matrix<double>::FromFile("data/h2o_sto-3g_two_electron_horton.data", rows, cols), std::runtime_error);  // can't read in two-electron data in a matrix


    // Test the read function on a small example mimicking the one-electron integrals
    GQCP::MatrixX<double> M_ref (rows, cols);
    M_ref << 2.1,  1.1,
             1.1, -3.4;

    auto M = GQCP::Matrix<double>::FromFile("data/small_one_ints.data", rows, cols);

    BOOST_CHECK(M.isApprox(M_ref, 1.0e-8));
}


BOOST_AUTO_TEST_CASE ( print ) {


    GQCP::MatrixX<double> M = GQCP::MatrixX<double>::Random(2, 2);

    std::ofstream file;
    file.open("print_output_stream_test.output");

    M.print();  // to std::cout
    M.print(file);

    file.close();
}


BOOST_AUTO_TEST_CASE ( minors ) {

    GQCP::MatrixX<double> A (3, 4);
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;


    GQCP::MatrixX<double> A_00 (2, 3);
    A_00 <<  6,  7,  8,
            10, 11, 12;
    BOOST_CHECK(A_00.isApprox(A.minor(0, 0)));

    Eigen::MatrixXd A_21 (2, 3);
    A_21 << 1, 3, 4,
            5, 7, 8;
    BOOST_CHECK(A_21.isApprox(A.minor(2, 1)));
}
