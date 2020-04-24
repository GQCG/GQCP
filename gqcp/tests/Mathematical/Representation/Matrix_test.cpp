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

#define BOOST_TEST_MODULE "Matrix"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Representation/Matrix.hpp"


BOOST_AUTO_TEST_CASE(constructor_assignment) {

    // A small check to see if the interface of the constructor and assignment operator works as expected

    GQCP::MatrixX<double> A = GQCP::MatrixX<double>::Random(3, 3);
    GQCP::MatrixX<double> B = GQCP::MatrixX<double>::Random(3, 3);

    GQCP::MatrixX<double> M1 {A * B};
    GQCP::MatrixX<double> M2 = A + B;
    GQCP::MatrixX<double> M3 = 2 * A;
}


BOOST_AUTO_TEST_CASE(Vector_FromFile) {

    size_t rows = 4;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.dat", rows), std::runtime_error);  // should be 'datA'


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.data", rows));


    // Check that there's an error when trying to read in tensor data into a vector
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("ref_data/h2o_sto-3g_two_electron_horton.data", rows), std::runtime_error);


    // Test the read function on a small example
    GQCP::VectorX<double> v_ref {rows};
    v_ref << 1.5, -0.2, 0.002, 8.3314;

    auto v = GQCP::VectorX<double>::FromFile("data/small_vector.data", rows);

    BOOST_CHECK(v.isApprox(v_ref, 1.0e-15));
}


BOOST_AUTO_TEST_CASE(Matrix_FromFile) {

    size_t rows = 2;
    size_t cols = 2;

    // Check that there's an error when a wrong path is supplied
    BOOST_CHECK_THROW(GQCP::Matrix<double>::FromFile("data/h2o_sto-3g_kinetic_horton.dat", rows, cols), std::runtime_error);  // should be 'datA'


    // Check that there's no error when a correct path is supplied
    BOOST_CHECK_NO_THROW(GQCP::Matrix<double>::FromFile("data/small_one_ints.data", rows, cols));


    // Check that there's an error when trying to read in tensor data into a matrix
    BOOST_CHECK_THROW(GQCP::Matrix<double>::FromFile("data/h2o_sto-3g_two_electron_horton.data", rows, cols), std::runtime_error);  // can't read in two-electron data in a matrix


    // Test the read function on a small example mimicking the one-electron integrals
    GQCP::MatrixX<double> M_ref {rows, cols};
    // clang-format off
    M_ref << 2.1,  1.1,
             1.1, -3.4;
    // clang-format on

    auto M = GQCP::Matrix<double>::FromFile("data/small_one_ints.data", rows, cols);

    BOOST_CHECK(M.isApprox(M_ref, 1.0e-8));
}


BOOST_AUTO_TEST_CASE(operator_call_CartesianDirection) {

    GQCP::Vector<size_t, 3> v = GQCP::Vector<size_t, 3>::Zero(3);
    v << 1, 2, 8;

    BOOST_CHECK(v(GQCP::CartesianDirection::x) == 1);
    BOOST_CHECK(v(GQCP::CartesianDirection::y) == 2);
    BOOST_CHECK(v(GQCP::CartesianDirection::z) == 8);
}


BOOST_AUTO_TEST_CASE(print) {

    GQCP::MatrixX<double> M = GQCP::MatrixX<double>::Random(2, 2);

    std::ofstream file;
    file.open("print_output_stream_test.output");

    M.print();  // to std::cout
    M.print(file);

    file.close();
}


BOOST_AUTO_TEST_CASE(matrixMinor) {

    GQCP::MatrixX<double> A {3, 4};
    // clang-format off
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;
    // clang-format on


    GQCP::MatrixX<double> A_00 {2, 3};
    // clang-format off
    A_00 <<  6,  7,  8,
            10, 11, 12;
    // clang-format on
    BOOST_CHECK(A_00.isApprox(A.matrixMinor(0, 0)));

    GQCP::MatrixX<double> A_21 {2, 3};
    // clang-format off
    A_21 << 1, 3, 4,
            5, 7, 8;
    // clang-format on
    BOOST_CHECK(A_21.isApprox(A.matrixMinor(2, 1)));
}


BOOST_AUTO_TEST_CASE(pairWiseReduce) {

    GQCP::MatrixX<double> M {4, 4};
    // clang-format off
    M <<  0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         12, 13, 14, 15;
    // clang-format on


    // Test default behavior
    GQCP::VectorX<double> v_ref_1 {16};
    v_ref_1 << 0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15;

    BOOST_CHECK(v_ref_1.isApprox(M.pairWiseReduce(), 1.0e-12));


    // Test non-default behavior
    GQCP::VectorX<double> v_ref_2 {9};
    v_ref_2 << 5, 9, 13, 6, 10, 14, 7, 11, 15;

    BOOST_CHECK(v_ref_2.isApprox(M.pairWiseReduce(1, 1), 1.0e-12));
}


/**
 *  Check the construction of a MatrixX from a given column-major vector
 */
BOOST_AUTO_TEST_CASE(FromColumnMajorVector) {

    // Construct a vector that is supposed to be in column-major ordering
    GQCP::VectorX<double> a {6};
    a << 1, 4, 2, 5, 3, 6;


    // Initialize the reference matrix and check the result
    GQCP::MatrixX<double> A_ref {2, 3};
    // clang-format off
    A_ref << 1, 2, 3,
             4, 5, 6;
    // clang-format on


    GQCP::MatrixX<double> A = GQCP::MatrixX<double>::FromColumnMajorVector(a, 2, 3);
    BOOST_CHECK(A.isApprox(A_ref, 1.0e-08));
}


/**
 *  Check the construction of a MatrixX from a given row-major vector
 */
BOOST_AUTO_TEST_CASE(FromRowMajorVector) {

    // Construct a vector that is supposed to be in column-major ordering
    GQCP::VectorX<double> a {6};
    a << 1, 2, 3, 4, 5, 6;


    // Initialize the reference matrix and check the result
    GQCP::MatrixX<double> A_ref {2, 3};
    // clang-format off
    A_ref << 1, 2, 3,
             4, 5, 6;
    // clang-format on


    GQCP::MatrixX<double> A = GQCP::MatrixX<double>::FromRowMajorVector(a, 2, 3);
    BOOST_CHECK(A.isApprox(A_ref, 1.0e-08));
}
