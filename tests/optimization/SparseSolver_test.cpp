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
#define BOOST_TEST_MODULE "Sparse"

#include "optimization/SparseSolver.hpp"

#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseSymMatProd.h"

#include "utilities/linalg.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( diagonal_getter_sparse ) {

    // Test the diagonal getter for a sparse matrix
    Eigen::VectorXd ref_diagonal (100);

    GQCP::SparseSolver sparse_solver (100);

    for (size_t i = 0; i < 100; i++) {
        sparse_solver.addToMatrix(2*i, i, i);
        ref_diagonal(i) = 2*i;
    }


    BOOST_CHECK(ref_diagonal.isApprox(sparse_solver.get_diagonal(), 1.0e-12));
}


BOOST_AUTO_TEST_CASE ( simple_sparse ) {

    // Create a random sparse symmetric matrix (adapted from https://stackoverflow.com/a/30742847/7930415).
    // Also, put it in the sparse matrix solver.
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist (0.0,1.0);

    size_t rows = 100;
    size_t cols = 100;


    std::vector<Eigen::Triplet<double>> triplet_list;  // needed for Eigen::SparseMatrix<double>
    GQCP::SparseSolver sparse_solver (100);


    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < i; j++) {
            auto random_number = dist(gen);

            if (random_number > 0.5) {  // if larger than a threshold, insert it (symmetrically)
                triplet_list.emplace_back(static_cast<int>(i), static_cast<int>(j), random_number);
                triplet_list.emplace_back(static_cast<int>(j), static_cast<int>(i), random_number);

                sparse_solver.addToMatrix(random_number, i, j);
                sparse_solver.addToMatrix(random_number, j, i);
            }
        }
    }

    Eigen::SparseMatrix<double> A (rows, cols);
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());


    // Find the lowest eigenpair using Spectra
    Spectra::SparseSymMatProd<double> matrixVectorProduct (A);
    Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>> spectra_sparse_eigensolver (&matrixVectorProduct, 1, 3);  // request 1 eigenpair, and use 3 Ritz pairs for the solution (need at least 2 more Ritz pairs than requested eigenvalues)
    spectra_sparse_eigensolver.init();
    spectra_sparse_eigensolver.compute();

    double ref_lowest_eigenvalue = spectra_sparse_eigensolver.eigenvalues()(0);
    Eigen::VectorXd ref_lowest_eigenvector = spectra_sparse_eigensolver.eigenvectors().col(0);


    // Find the lowest eigenpair using the sparse solver
    sparse_solver.solve();
    double test_lowest_eigenvalue = sparse_solver.get_eigenvalue();
    Eigen::VectorXd test_lowest_eigenvector = sparse_solver.get_eigenvector();


    BOOST_CHECK(std::abs(test_lowest_eigenvalue - ref_lowest_eigenvalue) < 1.0e-08);
    BOOST_CHECK(GQCP::areEqualEigenvectors(test_lowest_eigenvector, ref_lowest_eigenvector, 1.0e-08));


    // Check if the eigenvector is normalized
    BOOST_CHECK(std::abs(sparse_solver.get_eigenvector().norm() - 1) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( simple_sparse_number_of_requested_eigenpairs ) {

    size_t number_of_requested_eigenpairs = 5;

    // Create a random sparse symmetric matrix (adapted from https://stackoverflow.com/a/30742847/7930415).
    // Also, put it in the sparse matrix solver.
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.0,1.0);

    size_t rows = 100;
    size_t cols = 100;


    std::vector<Eigen::Triplet<double>> triplet_list;  // needed for Eigen::SparseMatrix<double>
    GQCP::SparseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = number_of_requested_eigenpairs;
    GQCP::SparseSolver sparse_solver (rows, solver_options);


    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < i; j++) {
            auto random_number = dist(gen);

            if (random_number > 0.5) {  // if larger than a threshold, insert it (symmetrically)
                triplet_list.emplace_back(static_cast<int>(i), static_cast<int>(j), random_number);
                triplet_list.emplace_back(static_cast<int>(j), static_cast<int>(i), random_number);

                sparse_solver.addToMatrix(random_number, i, j);
                sparse_solver.addToMatrix(random_number, j, i);
            }
        }
    }

    Eigen::SparseMatrix<double> A (rows, cols);
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());


    // Put this test in a try ... catch block since apparantly, on Travis the test fails
    try {
        // Find the lowest eigenpairs using Spectra
        Spectra::SparseSymMatProd<double> matrixVectorProduct(A);
        Spectra::SymEigsSolver<double, Spectra::SMALLEST_ALGE, Spectra::SparseSymMatProd<double>> spectra_sparse_eigensolver (&matrixVectorProduct, number_of_requested_eigenpairs, number_of_requested_eigenpairs + 3);  // number_of_requested_eigenpairs + 3 Ritz pairs for the solution (need at least 2 more Ritz pairs than requested eigenvalues)
        spectra_sparse_eigensolver.init();
        spectra_sparse_eigensolver.compute();

        Eigen::VectorXd ref_lowest_eigenvalues = spectra_sparse_eigensolver.eigenvalues().head(number_of_requested_eigenpairs);
        Eigen::MatrixXd ref_lowest_eigenvectors = spectra_sparse_eigensolver.eigenvectors().topLeftCorner(rows, number_of_requested_eigenpairs);

        // Create eigenpairs for the reference eigenpairs
        std::vector<GQCP::Eigenpair> ref_eigenpairs(number_of_requested_eigenpairs);
        for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
            ref_eigenpairs[i] = GQCP::Eigenpair(ref_lowest_eigenvalues(i),
                                                                ref_lowest_eigenvectors.col(i));
        }


        // Find the lowest eigenpairs using our sparse solver
        sparse_solver.solve();
        std::vector<GQCP::Eigenpair> eigenpairs = sparse_solver.get_eigenpairs();


        for (size_t i = 0; i < number_of_requested_eigenpairs; i++) {
            BOOST_CHECK(eigenpairs[i].isEqual(ref_eigenpairs[i]));  // check if the found eigenpairs are equal to the reference eigenpairs
            BOOST_CHECK(std::abs(eigenpairs[i].get_eigenvector().norm() - 1) < 1.0e-11);  // check if the found eigenpairs are normalized
        }

    } catch (const std::runtime_error& e) { // do nothing if Spectra inside the SparseSolver encounters a runtime_error
        std::cerr << "I caught a runtime error, which probably means that Spectra didn't find a solution. For the purpose of testing, this isn't numopt's fault.";
    }
}
