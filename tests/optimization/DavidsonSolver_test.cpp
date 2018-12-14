// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "DavidsonSolver"

#include "optimization/DavidsonSolver.hpp"

#include "utilities/linalg.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( constructor ) {

    // Create an example matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);

    // Test constructors with one supplied guess vector
    Eigen::VectorXd x_0 = Eigen::VectorXd::Constant(2, 1);

    // Create the solver options
    GQCP::DavidsonSolverOptions solver_options_faulty (x_0);
    solver_options_faulty.number_of_requested_eigenpairs = 3;
    solver_options_faulty.maximum_subspace_dimension = 4;
    solver_options_faulty.collapsed_subspace_dimension = 8;
    BOOST_CHECK_THROW(GQCP::DavidsonSolver (A, solver_options_faulty), std::invalid_argument);  // 3 requested eigenpairs: not enough initial guesses

    solver_options_faulty.number_of_requested_eigenpairs = 1;
    BOOST_CHECK_THROW(GQCP::DavidsonSolver (A, solver_options_faulty), std::invalid_argument);  // collapsed subspace dimension (8) cannot be larger than maximum subspace dimension (4)

    GQCP::DavidsonSolverOptions solver_options (x_0);
    BOOST_CHECK_NO_THROW(GQCP::DavidsonSolver (A, solver_options));


    // Test a constructor with two supplied guess vectors
    Eigen::MatrixXd Y_0 = Eigen::MatrixXd::Identity(2, 2);
    GQCP::DavidsonSolverOptions solver_options_faulty_2 (Y_0);
    solver_options_faulty_2.maximum_subspace_dimension = 8;
    solver_options_faulty_2.number_of_requested_eigenpairs = 2;
    solver_options_faulty_2.collapsed_subspace_dimension = 1;
    BOOST_CHECK_THROW(GQCP::DavidsonSolver (A, solver_options_faulty_2), std::invalid_argument);  // collapsed subspace dimension (1) cannot be smaller number of requested eigenpairs (2)
}


BOOST_AUTO_TEST_CASE ( constructor_raw ) {

    // Create an example matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);

    // Test constructors with one supplied guess vector
    Eigen::VectorXd x_0 = Eigen::VectorXd::Constant(2, 1);
    BOOST_CHECK_THROW(GQCP::DavidsonSolver davidson_solver (A, x_0, 3, 1.0e-08, 1.0e-12, 4, 8), std::invalid_argument);  // 3 requested eigenpairs: not enough initial guesses
    BOOST_CHECK_THROW(GQCP::DavidsonSolver davidson_solver (A, x_0, 1, 1.0e-08, 1.0e-12, 4, 8), std::invalid_argument);  // collapsed subspace dimension (8) cannot be larger than maximum subspace dimension (4)
    BOOST_CHECK_NO_THROW(GQCP::DavidsonSolver davidson_solver (A, x_0));


    // Test a constructor with two supplied guess vectors
    Eigen::MatrixXd Y_0 = Eigen::MatrixXd::Identity(2, 2);
    BOOST_CHECK_THROW(GQCP::DavidsonSolver davidson_solver (A, Y_0, 2, 1.0e-08, 1.0e-12, 8, 1), std::invalid_argument);  // collapsed subspace dimension (1) cannot be smaller number of requested eigenpairs (2)
}


BOOST_AUTO_TEST_CASE ( diagonal_getter_Davidson ) {

    // Test the diagonal getter for Davidson
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(2, 2);
    Eigen::VectorXd x_0 = Eigen::VectorXd::Constant(2, 1);

    GQCP::DavidsonSolver davidson_solver (A, GQCP::DavidsonSolverOptions (x_0));

    BOOST_CHECK(A.diagonal().isApprox(davidson_solver.get_diagonal(), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( esqc_example_solver ) {

    // We can find the following example at (http://www.esqc.org/static/lectures/Malmqvist_2B.pdf)


    // Build up the example matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Constant(5, 5, 0.1);
    A(0,0) = 1.0;
    A(1,1) = 2.0;
    A(2,2) = 3.0;
    A(3,3) = 3.0;
    A(4,4) = 3.0;


    // The solutions to the problem are given in the example
    double ref_lowest_eigenvalue = 0.979;

    Eigen::VectorXd ref_lowest_eigenvector (5);
    ref_lowest_eigenvector << 0.994, -0.083, -0.042, -0.042, -0.042;


    // Solve using the Davidson diagonalization, supplying an initial guess
    Eigen::VectorXd x_0 (5);
    x_0 << 1, 0, 0, 0, 0;
    GQCP::DavidsonSolver davidson_solver (A, GQCP::DavidsonSolverOptions (x_0));
    davidson_solver.solve();

    double test_lowest_eigenvalue = davidson_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = davidson_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 0.005);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_lowest_eigenvector, 0.005));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(davidson_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


// 12 iterations
BOOST_AUTO_TEST_CASE ( liu_50 ) {

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 50;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    double ref_lowest_eigenvalue = eigensolver.eigenvalues()(0);
    Eigen::VectorXd ref_lowest_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using the Davidson diagonalization, supplying an initial guess
    Eigen::VectorXd x_0 = Eigen::VectorXd::Zero(N);
    x_0(0) = 1;
    GQCP::DavidsonSolver davidson_solver (A, GQCP::DavidsonSolverOptions (x_0));
    davidson_solver.solve();

    double test_lowest_eigenvalue = davidson_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = davidson_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_lowest_eigenvector, 1.0e-08));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(davidson_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


// Test a forced subspace collapse for small dimensions (maximum subspace dimension < 12)
BOOST_AUTO_TEST_CASE ( liu_50_collapse ) {

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 50;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    double ref_eigenvalue = eigensolver.eigenvalues()(0);
    Eigen::VectorXd ref_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using the Davidson diagonalization, supplying an initial guess
    Eigen::VectorXd x_0 = Eigen::VectorXd::Zero(N);
    x_0(0) = 1;

    GQCP::DavidsonSolverOptions solver_options (x_0);
    solver_options.maximum_subspace_dimension = 10;
    solver_options.collapsed_subspace_dimension = 5;

    GQCP::DavidsonSolver davidson_solver (A, solver_options);  // maximum subspace dimension = 10, collapsed subspace dimension = 5
    davidson_solver.solve();

    double test_lowest_eigenvalue = davidson_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = davidson_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_eigenvalue) < 1.0e-08);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_eigenvector, 1.0e-08));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(davidson_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( liu_50_number_of_requested_eigenpairs ) {

    size_t number_of_requested_eigenpairs = 3;

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 50;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    Eigen::VectorXd ref_lowest_eigenvalues = eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    Eigen::MatrixXd ref_lowest_eigenvectors = eigensolver.eigenvectors().topLeftCorner(N, number_of_requested_eigenpairs);

    // Create eigenpairs for the reference eigenpairs
    std::vector<GQCP::Eigenpair> ref_eigenpairs (number_of_requested_eigenpairs);
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        ref_eigenpairs[i] = GQCP::Eigenpair(ref_lowest_eigenvalues(i), ref_lowest_eigenvectors.col(i));
    }

    // Solve using the Davidson diagonalization, supplying the requested amount of  initial guesses
    Eigen::MatrixXd X_0 = Eigen::MatrixXd::Identity(N, N).topLeftCorner(N, number_of_requested_eigenpairs);
    GQCP::DavidsonSolverOptions solver_options (X_0);
    solver_options.number_of_requested_eigenpairs = number_of_requested_eigenpairs;
    solver_options.collapsed_subspace_dimension = number_of_requested_eigenpairs;
    GQCP::DavidsonSolver davidson_solver (A, solver_options);
    davidson_solver.solve();

    std::vector<GQCP::Eigenpair> eigenpairs = davidson_solver.get_eigenpairs();



    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(eigenpairs[i].isEqual(ref_eigenpairs[i]));  // check if the found eigenpairs are equal to the reference eigenpairs
        BOOST_CHECK(std::abs(eigenpairs[i].get_eigenvector().norm() - 1) < 1.0e-12);  // check if the found eigenpairs are normalized
    }
}


// 12 iterations for large dimensions
BOOST_AUTO_TEST_CASE ( liu_1000 ) {

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 1000;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    double ref_eigenvalue = eigensolver.eigenvalues()(0);
    Eigen::VectorXd ref_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using the Davidson diagonalization, supplying an initial guess
    Eigen::VectorXd x_0 = Eigen::VectorXd::Zero(N);
    x_0(0) = 1;
    GQCP::DavidsonSolver davidson_solver (A, GQCP::DavidsonSolverOptions (x_0));
    davidson_solver.solve();

    double test_lowest_eigenvalue = davidson_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = davidson_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_eigenvalue) < 1.0e-08);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_eigenvector, 1.0e-08));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(davidson_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


// Test a forced subspace collapse for large dimensions (maximum subspace dimension < 12)
BOOST_AUTO_TEST_CASE ( liu_1000_collapse ) {

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 1000;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    double ref_eigenvalue = eigensolver.eigenvalues()(0);
    Eigen::VectorXd ref_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using the Davidson diagonalization, supplying an initial guess
    Eigen::VectorXd x_0 = Eigen::VectorXd::Zero(N);
    x_0(0) = 1;
    GQCP::DavidsonSolverOptions solver_options (x_0);
    solver_options.maximum_subspace_dimension = 10;
    solver_options.collapsed_subspace_dimension = 5;
    GQCP::DavidsonSolver davidson_solver (A, solver_options);  // maximum subspace dimension = 10, collapsed subspace dimension = 5
    davidson_solver.solve();

    double test_lowest_eigenvalue = davidson_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = davidson_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_eigenvalue) < 1.0e-08);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_eigenvector, 1.0e-08));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(davidson_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( liu_1000_number_of_requested_eigenpairs ) {

    size_t number_of_requested_eigenpairs = 3;

    // Let's prepare the Liu reference test (liu1978)
    size_t N = 1000;
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);
    Eigen::VectorXd ref_lowest_eigenvalues = eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    Eigen::MatrixXd ref_lowest_eigenvectors = eigensolver.eigenvectors().topLeftCorner(N, number_of_requested_eigenpairs);

    // Create eigenpairs for the reference eigenpairs
    std::vector<GQCP::Eigenpair> ref_eigenpairs (number_of_requested_eigenpairs);
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        ref_eigenpairs[i] = GQCP::Eigenpair(ref_lowest_eigenvalues(i), ref_lowest_eigenvectors.col(i));
    }

    // Solve using the Davidson diagonalization, supplying the requested amount of  initial guesses
    Eigen::MatrixXd X_0 = Eigen::MatrixXd::Identity(N, N).topLeftCorner(N, number_of_requested_eigenpairs);
    GQCP::DavidsonSolverOptions solver_options (X_0);
    solver_options.number_of_requested_eigenpairs = number_of_requested_eigenpairs;
    solver_options.collapsed_subspace_dimension = number_of_requested_eigenpairs;
    GQCP::DavidsonSolver davidson_solver (A, solver_options);
    davidson_solver.solve();

    std::vector<GQCP::Eigenpair> eigenpairs = davidson_solver.get_eigenpairs();



    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(eigenpairs[i].isEqual(ref_eigenpairs[i]));  // check if the found eigenpairs are equal to the reference eigenpairs
        BOOST_CHECK(std::abs(eigenpairs[i].get_eigenvector().norm() - 1) < 1.0e-12);  // check if the found eigenpairs are normalized
    }
}
