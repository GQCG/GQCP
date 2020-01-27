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
#define BOOST_TEST_MODULE "EigenproblemSolver"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"


/**
 *  Test if the eigenvalue problem solvers produce correct results.
 */


/**
 *  Test if the dense diagonalization algorithm finds the correct results. Even though it is a wrapper around Eigen's routines, we should check that the interplay between the corresponding environment is correctly implemented.
 */
BOOST_AUTO_TEST_CASE ( dense ) {

    const size_t number_of_requested_eigenpairs = 3;

    // Construct a random symmetric matrix
    const size_t dim = 10;
    GQCP::SquareMatrix<double> A = GQCP::SquareMatrix<double>::Random(dim,dim);
    GQCP::SquareMatrix<double> AT = A.transpose();
    GQCP::SquareMatrix<double> B = A + AT;


    // Diagonalize the matrix using Eigen's interface
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver (B);
    GQCP::VectorX<double> ref_lowest_eigenvalues = self_adjoint_eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
    GQCP::MatrixX<double> ref_lowest_eigenvectors = self_adjoint_eigensolver.eigenvectors().topLeftCorner(dim, number_of_requested_eigenpairs);

    // Create eigenpairs for the reference eigenpairs
    std::vector<GQCP::Eigenpair> ref_eigenpairs (number_of_requested_eigenpairs);
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        ref_eigenpairs[i] = GQCP::Eigenpair(ref_lowest_eigenvalues(i), ref_lowest_eigenvectors.col(i));
    }


    // Use our dense diagonalization algorithm to find the number of requested eigenpairs
    auto dense_environment = GQCP::EigenproblemEnvironment::Dense(B);
    auto dense_diagonalizer = GQCP::EigenproblemSolver::Dense(number_of_requested_eigenpairs);
    dense_diagonalizer.perform(dense_environment);

    const auto& eigenpairs = dense_environment.eigenpairs;


    // Check if the results are equal
    BOOST_CHECK(eigenpairs.size() == number_of_requested_eigenpairs);
    for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
        BOOST_CHECK(eigenpairs[i].isEqual(ref_eigenpairs[i]));  // check if the found eigenpairs are equal to the reference eigenpairs
        BOOST_CHECK(std::abs(eigenpairs[i].get_eigenvector().norm() - 1) < 1.0e-12);  // check if the found eigenpairs are normalized
    }
}
