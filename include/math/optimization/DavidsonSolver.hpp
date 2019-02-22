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
#ifndef GQCP_DAVIDSONSOLVER_HPP
#define GQCP_DAVIDSONSOLVER_HPP


#include "math/optimization/BaseEigenproblemSolver.hpp"
#include "math/optimization/EigenproblemSolverOptions.hpp"

#include "typedefs.hpp"



namespace GQCP {

/**
 *  A class that implements the Davidson algorithm for finding the lowest eigenpair of a (possibly large) diagonally-
 *  dominant symmetric matrix
 */
class DavidsonSolver : public BaseEigenproblemSolver {
private:
    double convergence_threshold;  // the tolerance on the norm of the residual vector
    double correction_threshold;  // the threshold used in solving the (approximated) residue correction equation
    size_t maximum_subspace_dimension;
    size_t collapsed_subspace_dimension;
    size_t maximum_number_of_iterations;
    size_t number_of_iterations = 0;

    VectorFunction matrixVectorProduct;
    Eigen::VectorXd diagonal;  // the diagonal of the matrix in question
    Eigen::MatrixXd V_0;  // the set of initial guesses (every column is an initial guess)


public:
    // CONSTRUCTORS
    /**
     *  @param matrixVectorProduct                  a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                             the diagonal of the matrix
     *  @param V_0                                  the (set of) initial guess(es) specified as a vector (matrix of column vectors)
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     *  @param convergence_threshold                the tolerance on the norm of the residual vector
     *  @param correction_threshold                 the threshold used in solving the (approximated) residue correction equation
     *  @param maximum_subspace_dimension           the maximum dimension of the Davidson subspace before collapsing
     *  @param collapsed_subspace_dimension         the dimension of the subspace after collapse
     *  @param maximum_number_of_iterations         the maximum number of Davidson iterations
     */
    DavidsonSolver(const VectorFunction& matrixVectorProduct, const Eigen::VectorXd& diagonal, const Eigen::MatrixXd& V_0, size_t number_of_requested_eigenpairs = 1, double convergence_threshold = 1.0e-08, double correction_threshold = 1.0e-12, size_t maximum_subspace_dimension = 15, size_t collapsed_subspace_dimension = 2, size_t maximum_number_of_iterations = 128);

    /**
     *  @param A                                    the matrix to be diagonalized
     *  @param V_0                                  the (set of) initial guess(es) specified as a vector (matrix of column vectors)
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     *  @param convergence_threshold                the tolerance on the norm of the residual vector
     *  @param correction_threshold                 the threshold used in solving the (approximated) residue correction equation
     *  @param maximum_subspace_dimension           the maximum dimension of the Davidson subspace before collapsing
     *  @param collapsed_subspace_dimension         the dimension of the subspace after collapse
     *  @param maximum_number_of_iterations         the maximum number of Davidson iterations
     */
    DavidsonSolver(const Eigen::MatrixXd& A, const Eigen::MatrixXd& V_0, size_t number_of_requested_eigenpairs = 1, double convergence_threshold = 1.0e-08, double correction_threshold = 1.0e-12, size_t maximum_subspace_dimension = 15, size_t collapsed_subspace_dimension = 2, size_t maximum_number_of_iterations = 128);

    /**
     *  @param matrixVectorProduct          a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                     the diagonal of the matrix
     *  @param davidson_solver_options      the options specified for solving the Davidson eigenvalue problem
     */
    DavidsonSolver(const VectorFunction& matrixVectorProduct, const Eigen::VectorXd& diagonal, const DavidsonSolverOptions& davidson_solver_options);

    /**
     *  @param A                            the matrix to be diagonalized
     *  @param davidson_solver_options      the options specified for solving the Davidson eigenvalue problem
     */
    DavidsonSolver(const Eigen::MatrixXd& A, const DavidsonSolverOptions& davidson_solver_options);


    // DESTRUCTOR
    ~DavidsonSolver() override = default;


    // GETTERS
    const Eigen::VectorXd& get_diagonal() const { return this->diagonal; };
    size_t get_number_of_iterations() const;


    // PUBLIC METHODS
    /**
     *  Solve the eigenvalue problem related to the given matrix-vector product
     *
     *  If successful, it sets
     *      - _is_solved to true
     *      - the number of requested eigenpairs
     */
    void solve() override;
};


}  // namespace GQCP



#endif  // GQCP_DAVIDSONSOLVER_HPP
