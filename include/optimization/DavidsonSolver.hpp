// This file is part of GQCG-numopt.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-numopt is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-numopt is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-numopt.  If not, see <http://www.gnu.org/licenses/>.
#ifndef NUMOPT_DAVIDSONSOLVER_HPP
#define NUMOPT_DAVIDSONSOLVER_HPP


#include "BaseEigenproblemSolver.hpp"
#include "EigenproblemSolverOptions.hpp"

#include "common.hpp"



namespace numopt {
namespace eigenproblem {

/**
 *  A class that implements the Davidson algorithm for finding the lowest eigenpair of a (possibly large) diagonally-
 *  dominant symmetric matrix
 */
class DavidsonSolver : public numopt::eigenproblem::BaseEigenproblemSolver {
private:
    static constexpr size_t maximum_number_of_iterations = 128;

    const double convergence_threshold;  // the tolerance on the norm of the residual vector
    const double correction_threshold;  // the threshold used in solving the (approximated) residue correction equation
    const size_t maximum_subspace_dimension;
    const size_t collapsed_subspace_dimension;

    const numopt::VectorFunction matrixVectorProduct;
    const Eigen::VectorXd diagonal;  // the diagonal of the matrix in question
    const Eigen::MatrixXd V_0;  // the set of initial guesses (every column is an initial guess)

    size_t number_of_iterations = 0;



public:
    // CONSTRUCTORS
    /**
     *  @param matrixVectorProduct          a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                     the diagonal of the matrix
     *  @param davidson_solver_options      the options specified for solving the Davidson eigenvalue problem
     */
    DavidsonSolver(const numopt::VectorFunction& matrixVectorProduct, const Eigen::VectorXd& diagonal, const DavidsonSolverOptions& davidson_solver_options);

    /**
     *  @param A                            the matrix to be diagonalized
     *  @param davidson_solver_options      the options specified for solving the Davidson eigenvalue problem
     */
    DavidsonSolver(const Eigen::MatrixXd& A, const DavidsonSolverOptions& davidson_solver_options);

    /**
     *  @param matrixVectorProduct                  a vector function that returns the matrix-vector product (i.e. the matrix-vector product representation of the matrix)
     *  @param diagonal                             the diagonal of the matrix
     *  @param V_0                                  the (set of) initial guess(es) specified as a vector (matrix of column vectors)
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     *  @param residue_tolerance                    the tolerance on the norm of the residual vector
     *  @param correction_threshold                 the threshold used in solving the (approximated) residue correction equation
     *  @param maximum_subspace_dimension           the maximum dimension of the Davidson subspace before collapsing
     *  @param collapsed_subspace_dimension         the dimension of the subspace after collapse
     */
    DavidsonSolver(const numopt::VectorFunction& matrixVectorProduct, const Eigen::VectorXd& diagonal,
                   const Eigen::MatrixXd& V_0, size_t number_of_requested_eigenpairs = 1,
                   double residue_tolerance = 1.0e-08, double correction_threshold = 1.0e-12,
                   size_t maximum_subspace_dimension = 15, size_t collapsed_subspace_dimension = 2);

    /**
     *  @param A                                    the matrix to be diagonalized
     *  @param V_0                                  the (set of) initial guess(es) specified as a vector (matrix of column vectors)
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the solver should find
     *  @param residue_tolerance                    the tolerance on the norm of the residual vector
     *  @param correction_threshold                 the threshold used in solving the (approximated) residue correction equation
     *  @param maximum_subspace_dimension           the maximum dimension of the Davidson subspace before collapsing
     *  @param collapsed_subspace_dimension         the dimension of the subspace after collapse
     */
    DavidsonSolver(const Eigen::MatrixXd& A, const Eigen::MatrixXd& V_0, size_t number_of_requested_eigenpairs = 1,
                   double residue_tolerance = 1.0e-08, double correction_threshold = 1.0e-12,
                   size_t maximum_subspace_dimension = 15, size_t collapsed_subspace_dimension = 2);


    // DESTRUCTOR
    ~DavidsonSolver() override = default;


    // GETTERS
    Eigen::VectorXd get_diagonal() override { return this->diagonal; };
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


}  // namespace eigenproblem
}  // namespace numopt



#endif  // NUMOPT_DAVIDSONSOLVER_HPP
