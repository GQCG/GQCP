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

#define BOOST_TEST_MODULE "EigenproblemSolver"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Mathematical/Representation/Matrix.hpp"


/**
 *  Test if the eigenvalue problem solvers produce correct results.
 */


/**
 *  Test if the dense diagonalization algorithm finds the correct results. Even though it is a wrapper around Eigen's routines, we should check that the interplay between the solver and its corresponding environment is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(dense) {

    const size_t number_of_requested_eigenpairs = 3;

    // Construct a random symmetric matrix.
    const size_t dim = 10;
    const GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Random(dim);
    const GQCP::SquareMatrix<double> AT = A.transpose();
    const GQCP::SquareMatrix<double> B = A + AT;


    // Diagonalize the matrix using Eigen's interface
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver {B};
    const GQCP::VectorX<double> ref_lowest_eigenvalues = self_adjoint_eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    const GQCP::MatrixX<double> ref_lowest_eigenvectors = self_adjoint_eigensolver.eigenvectors().topLeftCorner(dim, number_of_requested_eigenpairs);

    // Create eigenpairs for the reference eigenpairs.
    std::vector<GQCP::Eigenpair<double>> ref_eigenpairs {number_of_requested_eigenpairs};
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        ref_eigenpairs[i] = GQCP::Eigenpair<double>(ref_lowest_eigenvalues(i), ref_lowest_eigenvectors.col(i));
    }


    // Use our dense diagonalization algorithm to find the number of requested eigenpairs.
    auto dense_environment = GQCP::EigenproblemEnvironment::Dense(B);
    auto dense_diagonalizer = GQCP::EigenproblemSolver::Dense();
    dense_diagonalizer.perform(dense_environment);

    const auto eigenpairs = dense_environment.eigenpairs(number_of_requested_eigenpairs);


    // Check if the results are equal.
    BOOST_CHECK(eigenpairs.size() == number_of_requested_eigenpairs);
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(eigenpairs[i].isEqualTo(ref_eigenpairs[i]));                  // Check if the found eigenpairs are equal to the reference eigenpairs.
        BOOST_CHECK(std::abs(eigenpairs[i].eigenvector().norm() - 1) < 1.0e-12);  // Check if the found eigenpairs are normalized.
    }
}


/**
 *  Check if the Davidson algorithm finds correct results. The example is taken from (http://www.esqc.org/static/lectures/Malmqvist_2B.pdf).
 */
BOOST_AUTO_TEST_CASE(Davidson_ESQC_example_solver) {

    // Build up the example matrix.
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Constant(5, 5, 0.1);
    A(0, 0) = 1.0;
    A(1, 1) = 2.0;
    A(2, 2) = 3.0;
    A(3, 3) = 3.0;
    A(4, 4) = 3.0;


    // The solutions to the problem are given in the example
    const double ref_lowest_eigenvalue = 0.979;

    GQCP::VectorX<double> ref_lowest_eigenvector {5};
    ref_lowest_eigenvector << 0.994, -0.083, -0.042, -0.042, -0.042;


    // Solve using the Davidson diagonalization, supplying an initial guess.
    GQCP::VectorX<double> x_0 {5};
    x_0 << 1, 0, 0, 0, 0;
    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, x_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson();  // The default is finding only the eigenpair with the lowest eigenvalue.
    davidson_solver.perform(davidson_environment);

    const double test_lowest_eigenvalue = davidson_environment.eigenvalues(0);
    const GQCP::VectorX<double> test_lowest_eigenvector = davidson_environment.eigenvectors.col(0);


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 0.005);
    BOOST_CHECK(test_lowest_eigenvector.isEqualEigenvectorAs(ref_lowest_eigenvector, 0.005));
    BOOST_CHECK(std::abs(test_lowest_eigenvector.norm() - 1) < 1.0e-12);
}


/**
 *  Check the workings of the Davidson algorithm for Liu's reference test (article: Liu1978).
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_50) {

    // Build up the example matrix.
    const size_t N = 50;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const double ref_lowest_eigenvalue = eigensolver.eigenvalues()(0);
    const GQCP::VectorX<double> ref_lowest_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using our Davidson diagonalization algorithm, supplying an initial guess.
    GQCP::VectorX<double> x_0 = GQCP::VectorX<double>::Unit(N, 0);

    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, x_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson();  // The default is finding only the eigenpair with the lowest eigenvalue.
    davidson_solver.perform(davidson_environment);

    const double test_lowest_eigenvalue = davidson_environment.eigenvalues(0);
    const GQCP::VectorX<double> test_lowest_eigenvector = davidson_environment.eigenvectors.col(0);


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(test_lowest_eigenvector.isEqualEigenvectorAs(ref_lowest_eigenvector, 1.0e-08));
    BOOST_CHECK(std::abs(test_lowest_eigenvector.norm() - 1) < 1.0e-12);
}


/**
 *  Check if the Davidson algorithm works for Liu's reference test (Liu1975) when a subspace collapse is forced.
 * 
 *  This test uses the same example as the previous test, from which we know that it requires 12 iterations. In this test, we set the maximum subspace dimensions to 10, in order to force the subspace collapse.
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_50_collapse) {

    // Build up the example matrix.
    const size_t N = 50;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const double ref_lowest_eigenvalue = eigensolver.eigenvalues()(0);
    const GQCP::VectorX<double> ref_lowest_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using our Davidson diagonalization algorithm, supplying an initial guess
    GQCP::VectorX<double> x_0 = GQCP::VectorX<double>::Unit(N, 0);

    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, x_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson(1, 10);  // number_of_requested_eigenpairs=1, maximum_subspace_dimension=10
    davidson_solver.perform(davidson_environment);


    const double test_lowest_eigenvalue = davidson_environment.eigenvalues(0);
    const GQCP::VectorX<double> test_lowest_eigenvector = davidson_environment.eigenvectors.col(0);


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(test_lowest_eigenvector.isEqualEigenvectorAs(ref_lowest_eigenvector, 1.0e-08));
    BOOST_CHECK(std::abs(test_lowest_eigenvector.norm() - 1) < 1.0e-12);
}


/**
 *  Check if the Davidson algorithm works for Liu's reference test (Liu1978), when the number of requested eigenpairs is different from 1.
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_50_number_of_requested_eigenpairs) {

    const size_t number_of_requested_eigenpairs = 3;

    // Build up the example matrix.
    const size_t N = 50;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const GQCP::VectorX<double> ref_lowest_eigenvalues = eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    const GQCP::MatrixX<double> ref_lowest_eigenvectors = eigensolver.eigenvectors().topLeftCorner(N, number_of_requested_eigenpairs);


    // Solve using our Davidson diagonalization algorithm, supplying a number of initial guesses.
    const GQCP::MatrixX<double> X_0 = GQCP::MatrixX<double>::Identity(N, N).topLeftCorner(N, number_of_requested_eigenpairs);

    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, X_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson(3);  // The number of requested eigenpairs is 3.
    davidson_solver.perform(davidson_environment);


    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(std::abs(davidson_environment.eigenvalues(i) - ref_lowest_eigenvalues(i)));

        const GQCP::VectorX<double> davidson_eigenvector = davidson_environment.eigenvectors.col(i);
        const GQCP::VectorX<double> ref_eigenvector = ref_lowest_eigenvectors.col(i);
        BOOST_CHECK(davidson_eigenvector.isEqualEigenvectorAs(ref_eigenvector, 1.0e-08));

        BOOST_CHECK(std::abs(davidson_environment.eigenvectors.col(i).norm() - 1) < 1.0e-12);  // Check if the found eigenpairs are normalized.
    }
}


/**
 *  Check the workings of the Davidson algorithm for Liu's reference test (article: Liu1978) with large dimensions.
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_1000) {

    // Build up the example matrix
    const size_t N = 1000;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const double ref_lowest_eigenvalue = eigensolver.eigenvalues()(0);
    const GQCP::VectorX<double> ref_lowest_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using our Davidson diagonalization algorithm, supplying an initial guess.
    GQCP::VectorX<double> x_0 = GQCP::VectorX<double>::Unit(N, 0);
    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, x_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson();  // The default is finding only the eigenpair with the lowest eigenvalue.
    davidson_solver.perform(davidson_environment);

    const double test_lowest_eigenvalue = davidson_environment.eigenvalues(0);
    const GQCP::VectorX<double> test_lowest_eigenvector = davidson_environment.eigenvectors.col(0);


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(test_lowest_eigenvector.isEqualEigenvectorAs(ref_lowest_eigenvector, 1.0e-08));
    BOOST_CHECK(std::abs(test_lowest_eigenvector.norm() - 1) < 1.0e-12);
}


/**
 *  Check if the Davidson algorithm works when a subspace collapse is forced with large dimensions.
 * 
 *  This test uses the same example as the previous test, from which we know that it requires 12 iterations. In this test, we set the maximum subspace dimensions to 10, in order to force the subspace collapse.
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_1000_collapse) {

    // Build up the example matrix.
    const size_t N = 1000;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const double ref_lowest_eigenvalue = eigensolver.eigenvalues()(0);
    const GQCP::VectorX<double> ref_lowest_eigenvector = eigensolver.eigenvectors().col(0);


    // Solve using our Davidson diagonalization algorithm, supplying an initial guess
    GQCP::VectorX<double> x_0 = GQCP::VectorX<double>::Unit(N, 0);

    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, x_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson(1, 10);  // number_of_requested_eigenpairs=1, maximum_subspace_dimension=10
    davidson_solver.perform(davidson_environment);


    const double test_lowest_eigenvalue = davidson_environment.eigenvalues(0);
    const GQCP::VectorX<double> test_lowest_eigenvector = davidson_environment.eigenvectors.col(0);


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(test_lowest_eigenvector.isEqualEigenvectorAs(ref_lowest_eigenvector, 1.0e-08));
    BOOST_CHECK(std::abs(test_lowest_eigenvector.norm() - 1) < 1.0e-12);
}


/**
 *  Check if the Davidson algorithm works for a number of requested eigenpairs different from 1 with large dimensions.
 */
BOOST_AUTO_TEST_CASE(Davidson_Liu_1000_number_of_requested_eigenpairs) {

    const size_t number_of_requested_eigenpairs = 3;

    // Build up the example matrix.
    const size_t N = 1000;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Ones(N, N);
    for (size_t i = 0; i < N; i++) {
        if (i < 5) {
            A(i, i) = 1 + 0.1 * i;
        } else {
            A(i, i) = 2 * (i + 1) - 1;
        }
    }


    // Solve the eigenvalue problem with Eigen.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {A};
    const GQCP::VectorX<double> ref_lowest_eigenvalues = eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    const GQCP::MatrixX<double> ref_lowest_eigenvectors = eigensolver.eigenvectors().topLeftCorner(N, number_of_requested_eigenpairs);


    // Solve using our Davidson diagonalization algorithm, supplying a number of initial guesses.
    const GQCP::MatrixX<double> X_0 = GQCP::MatrixX<double>::Identity(N, N).topLeftCorner(N, number_of_requested_eigenpairs);

    auto davidson_environment = GQCP::EigenproblemEnvironment::Iterative(A, X_0);
    auto davidson_solver = GQCP::EigenproblemSolver::Davidson(3, 10);  // number_of_requested_eigenpairs=3, maximum_subspace_dimension=10
    davidson_solver.perform(davidson_environment);


    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(std::abs(davidson_environment.eigenvalues(i) - ref_lowest_eigenvalues(i)));

        const GQCP::VectorX<double> davidson_eigenvector = davidson_environment.eigenvectors.col(i);
        const GQCP::VectorX<double> ref_eigenvector = ref_lowest_eigenvectors.col(i);
        BOOST_CHECK(davidson_eigenvector.isEqualEigenvectorAs(ref_eigenvector, 1.0e-08));

        BOOST_CHECK(std::abs(davidson_environment.eigenvectors.col(i).norm() - 1) < 1.0e-12);
    }
}
