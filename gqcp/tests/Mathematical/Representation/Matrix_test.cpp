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


/**
 *  Test the named constructor `FromFile` for VectorX.
 */
BOOST_AUTO_TEST_CASE(Vector_FromFile) {

    const size_t rows = 4;

    // Check that there's an error when a wrong path is supplied.
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.dat", rows), std::runtime_error);  // Should be 'datA'.


    // Check that there's no error when a correct path is supplied.
    BOOST_CHECK_NO_THROW(GQCP::VectorX<double>::FromFile("data/small_vector.data", rows));


    // Check that there's an error when trying to read in tensor data into a vector.
    BOOST_CHECK_THROW(GQCP::VectorX<double>::FromFile("ref_data/h2o_sto-3g_two_electron_horton.data", rows), std::runtime_error);


    // Test the read function on a small example.
    GQCP::VectorX<double> v_ref {rows};
    v_ref << 1.5, -0.2, 0.002, 8.3314;

    const auto v = GQCP::VectorX<double>::FromFile("data/small_vector.data", rows);

    BOOST_CHECK(v.isApprox(v_ref, 1.0e-15));
}


/**
 *  Test the named constructor `FromFile` for MatrixX.
 */
BOOST_AUTO_TEST_CASE(Matrix_FromFile) {

    size_t rows = 2;
    size_t cols = 2;

    // Check that there's an error when a wrong path is supplied.
    BOOST_CHECK_THROW(GQCP::MatrixX<double>::FromFile("data/h2o_sto-3g_kinetic_horton.dat", rows, cols), std::runtime_error);  // Should be 'datA'.


    // Check that there's no error when a correct path is supplied.
    BOOST_CHECK_NO_THROW(GQCP::MatrixX<double>::FromFile("data/small_one_ints.data", rows, cols));


    // Check that there's an error when trying to read in tensor data into a matrix.
    BOOST_CHECK_THROW(GQCP::MatrixX<double>::FromFile("data/h2o_sto-3g_two_electron_horton.data", rows, cols), std::runtime_error);  // Can't read in two-electron data in a matrix.


    // Test the read function on a small example mimicking the one-electron integrals.
    GQCP::MatrixX<double> M_ref {rows, cols};
    // clang-format off
    M_ref << 2.1,  1.1,
             1.1, -3.4;
    // clang-format on

    const auto M = GQCP::MatrixX<double>::FromFile("data/small_one_ints.data", rows, cols);

    BOOST_CHECK(M.isApprox(M_ref, 1.0e-8));
}


/**
 *  Check a Vector's `isEqualEigenvectorAs` method.
 */
BOOST_AUTO_TEST_CASE(isEqualEigenvectorAs) {

    // Test isEqualEigenvectorAs with an example.
    const GQCP::Vector<double, 3> a {2, 3, 1};
    const GQCP::Vector<double, 3> b {2, 3, 1};
    const GQCP::Vector<double, 3> c {-2, -3, -1};
    const GQCP::Vector<double, 3> d {2, 3, 0};

    BOOST_CHECK(a.isEqualEigenvectorAs(b, 1.0e-6));
    BOOST_CHECK(a.isEqualEigenvectorAs(c, 1.0e-6));  // Eigenvectors are equal up to their sign.
    BOOST_CHECK(b.isEqualEigenvectorAs(c, 1.0e-6));  // Eigenvectors are equal up to their sign.

    BOOST_CHECK(!a.isEqualEigenvectorAs(d, 1.0e-6));
    BOOST_CHECK(!c.isEqualEigenvectorAs(d, 1.0e-6));
}


/**
 *  Check the argument parsing for `hasEqualSetsOfEigenvectorsAs`.
 */
BOOST_AUTO_TEST_CASE(hasEqualSetsOfEigenvectorsAs_throws) {

    // Check for throws if the dimensions aren't compatible.
    const GQCP::MatrixX<double> C1 {3, 3};
    const GQCP::MatrixX<double> C2 {3, 2};

    BOOST_CHECK_THROW(C1.hasEqualSetsOfEigenvectorsAs(C2, 1.0e-6), std::invalid_argument);


    // Check for no throw if the dimensions are compatible
    const GQCP::MatrixX<double> C3 {3, 3};

    BOOST_CHECK_NO_THROW(C1.hasEqualSetsOfEigenvectorsAs(C3, 1.0e-6));
}


/**
 *  Check a Matrix' `hasEqualSetsOfEigenvectorsAs` method.
 */
BOOST_AUTO_TEST_CASE(hasEqualSetsOfEigenvectorsAs_example) {

    // Test hasEqualSetsOfEigenvectorsAs with an example.
    GQCP::MatrixX<double> eigenvectors1 {2, 2};
    GQCP::MatrixX<double> eigenvectors2 {2, 2};
    GQCP::MatrixX<double> eigenvectors3 {2, 2};
    GQCP::MatrixX<double> eigenvectors4 {2, 2};

    // clang-format off
    eigenvectors1 << 0,  2,
                     1, -1;

    eigenvectors2 << 0,  2,
                     1, -1;

    eigenvectors3 << 0, -2,
                     1,  1;

    eigenvectors4 << 1,  2,
                     1, -1;
    // clang-format on

    BOOST_CHECK(eigenvectors1.hasEqualSetsOfEigenvectorsAs(eigenvectors2, 1.0e-6));
    BOOST_CHECK(eigenvectors1.hasEqualSetsOfEigenvectorsAs(eigenvectors3, 1.0e-6));

    BOOST_CHECK(!eigenvectors1.hasEqualSetsOfEigenvectorsAs(eigenvectors4, 1.0e-6));
}


/**
 *  Check if `CartesianDirection` may be given to the call operator.
 */
BOOST_AUTO_TEST_CASE(operator_call_CartesianDirection) {

    const GQCP::Vector<size_t, 3> v {1, 2, 8};

    BOOST_CHECK(v(GQCP::CartesianDirection::x) == 1);
    BOOST_CHECK(v(GQCP::CartesianDirection::y) == 2);
    BOOST_CHECK(v(GQCP::CartesianDirection::z) == 8);
}


/**
 *  Check if a `MatrixX` may be printed.
 */
BOOST_AUTO_TEST_CASE(print) {

    const GQCP::MatrixX<double> M = GQCP::MatrixX<double>::Random(2, 2);

    std::ofstream file;
    file.open("print_output_stream_test.output");

    M.print();      // To std::cout.
    M.print(file);  // To the given file.

    file.close();
}


/**
 *  Check if the minor of a matrix is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(calculateMinor) {

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
    BOOST_CHECK(A_00.isApprox(A.calculateMinor(0, 0)));

    GQCP::MatrixX<double> A_21 {2, 3};
    // clang-format off
    A_21 << 1, 3, 4,
            5, 7, 8;
    // clang-format on
    BOOST_CHECK(A_21.isApprox(A.calculateMinor(2, 1)));
}


/**
 *  Check if pair-wise reshaping a matrix to a vector is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(pairWiseReduced) {

    GQCP::MatrixX<double> M {4, 4};
    // clang-format off
    M <<  0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         12, 13, 14, 15;
    // clang-format on


    // Test default behavior.
    GQCP::VectorX<double> v_ref_1 {16};
    v_ref_1 << 0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15;

    BOOST_CHECK(v_ref_1.isApprox(M.pairWiseReduced(), 1.0e-12));


    // Test non-default behavior.
    GQCP::VectorX<double> v_ref_2 {9};
    v_ref_2 << 5, 9, 13, 6, 10, 14, 7, 11, 15;

    BOOST_CHECK(v_ref_2.isApprox(M.pairWiseReduced(1, 1), 1.0e-12));
}


/**
 *  Check the construction of a MatrixX from a given column-major vector.
 */
BOOST_AUTO_TEST_CASE(FromColumnMajorVector) {

    // Construct a vector that is supposed to be in column-major ordering.
    GQCP::VectorX<double> a {6};
    a << 1, 4, 2, 5, 3, 6;


    // Initialize the reference matrix and check the result.
    GQCP::MatrixX<double> A_ref {2, 3};
    // clang-format off
    A_ref << 1, 2, 3,
             4, 5, 6;
    // clang-format on


    const auto A = GQCP::MatrixX<double>::FromColumnMajorVector(a, 2, 3);
    BOOST_CHECK(A.isApprox(A_ref, 1.0e-08));
}


/**
 *  Check the construction of a MatrixX from a given row-major vector.
 */
BOOST_AUTO_TEST_CASE(FromRowMajorVector) {

    // Construct a vector that is supposed to be in column-major ordering.
    GQCP::VectorX<double> a {6};
    a << 1, 2, 3, 4, 5, 6;


    // Initialize the reference matrix and check the result.
    GQCP::MatrixX<double> A_ref {2, 3};
    // clang-format off
    A_ref << 1, 2, 3,
             4, 5, 6;
    // clang-format on


    GQCP::MatrixX<double> A = GQCP::MatrixX<double>::FromRowMajorVector(a, 2, 3);
    BOOST_CHECK(A.isApprox(A_ref, 1.0e-08));
}


/**
 *  Check if removing single rows and columns works as expected.
 */
BOOST_AUTO_TEST_CASE(remove_single_rows_columns) {

    GQCP::MatrixX<double> M {4, 4};
    // clang-format off
    M <<  0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         12, 13, 14, 15;
    // clang-format on


    // Remove the row with index 2.
    GQCP::MatrixX<double> M_ref1 {3, 4};
    // clang-format off
    M_ref1 <<  0,  1,  2,  3,
               4,  5,  6,  7,
              12, 13, 14, 15;
    // clang-format on
    GQCP::MatrixX<double> M_1 = M;
    M_1.removeRow(2);
    BOOST_CHECK(M_1.isApprox(M_ref1, 1.0e-08));


    // Remove the column with index 1.
    GQCP::MatrixX<double> M_ref2 {4, 3};
    // clang-format off
    M_ref2 <<  0,  2,  3,
               4,  6,  7,
               8, 10, 11,
              12, 14, 15;
    // clang-format on
    GQCP::MatrixX<double> M_2 = M;
    M_2.removeColumn(1);
    BOOST_CHECK(M_2.isApprox(M_ref2, 1.0e-08));


    // Remove the column with index 3.
    GQCP::MatrixX<double> M_ref3 {4, 3};
    // clang-format off
    M_ref3 <<  0,  1,  2,
               4,  5,  6,
               8,  9, 10,
              12, 13, 14;
    // clang-format on
    GQCP::MatrixX<double> M_3 = M;
    M_3.removeColumn(3);
    BOOST_CHECK(M_3.isApprox(M_ref3, 1.0e-08));
}


/**
 *  Check if removing multiple rows and columns works as expected.
 */
BOOST_AUTO_TEST_CASE(remove_multiple_rows_columns) {

    GQCP::MatrixX<double> M {4, 4};
    // clang-format off
    M <<  0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         12, 13, 14, 15;
    // clang-format on


    // Remove the rows with indices 2 and 3.
    GQCP::MatrixX<double> M_ref1 {2, 4};
    // clang-format off
    M_ref1 <<  0, 1, 2, 3,
               4, 5, 6, 7;
    // clang-format on
    GQCP::MatrixX<double> M_1 = M;
    M_1.removeRows({2, 3});
    BOOST_CHECK(M_1.isApprox(M_ref1, 1.0e-08));


    // Remove the columns with indices 1 and 3.
    GQCP::MatrixX<double> M_ref2 {4, 2};
    // clang-format off
    M_ref2 <<  0,  2,
               4,  6,
               8, 10,
              12, 14;
    // clang-format on
    GQCP::MatrixX<double> M_2 = M;
    M_2.removeColumns({1, 3});
    BOOST_CHECK(M_2.isApprox(M_ref2, 1.0e-08));


    // Remove the column with index 3.
    GQCP::MatrixX<double> M_ref3 {4, 3};
    // clang-format off
    M_ref3 <<  0,  1,  2,
               4,  5,  6,
               8,  9, 10,
              12, 13, 14;
    // clang-format on
    GQCP::MatrixX<double> M_3 = M;
    M_3.removeColumns({3});
    BOOST_CHECK(M_3.isApprox(M_ref3, 1.0e-08));
}
