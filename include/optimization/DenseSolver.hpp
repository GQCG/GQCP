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
#ifndef NUMOPT_DENSESOLVER_HPP
#define NUMOPT_DENSESOLVER_HPP



#include "BaseMatrixSolver.hpp"
#include "EigenproblemSolverOptions.hpp"




namespace numopt {
namespace eigenproblem {


/**
 *  An eigenproblem solver that stores a dense representation of the matrix
 */
class DenseSolver : public numopt::eigenproblem::BaseMatrixSolver {
private:
    Eigen::MatrixXd matrix;


public:
    // CONSTRUCTORS
    /**
     *  @param dim                      the dimension of the matrix
     *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
     */
    DenseSolver(size_t dim, const DenseSolverOptions& dense_solver_options);

    /**
     *  @param matrix                   the full dense representation of the matrix
     *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
     */
    DenseSolver(const Eigen::MatrixXd& matrix, const DenseSolverOptions& dense_solver_options);

    /**
     *  @param dim                                  the dimension of the matrix
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
     */
    DenseSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);

    /**
     *  @param matrix                   the full dense representation of the matrix
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
     */
    DenseSolver(Eigen::MatrixXd matrix, size_t number_of_requested_eigenpairs = 1);


    // DESTRUCTOR
    ~DenseSolver() override = default;


    // GETTERS
    Eigen::MatrixXd get_matrix() { return this->matrix; };
    Eigen::VectorXd get_diagonal() override { return this->matrix.diagonal(); };


    // PUBLIC OVERRIDDEN METHODS
    /**
     *  Solve the full dense eigenvalue problem
     *
     *  If successful, it sets
     *      - _is_solved to true
     *      - the number of requested eigenpairs
     */
    void solve() override;

    /**
     *  @param value        the value to be added
     *  @param index1       the first index of the matrix
     *  @param index2       the second index of the matrix
     *
     *  Add the value to the matrix at (index1, index2)
     */
    void addToMatrix(double value, size_t index1, size_t index2) override;
};



}  // namespace eigenproblem
}  // namespace numopt



#endif  // NUMOPT_DENSESOLVER_HPP
