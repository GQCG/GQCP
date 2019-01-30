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
#ifndef GQCP_DENSESOLVER_HPP
#define GQCP_DENSESOLVER_HPP



#include "optimization/BaseMatrixSolver.hpp"
#include "optimization/EigenproblemSolverOptions.hpp"




namespace GQCP {


/**
 *  An eigenproblem solver that stores a dense representation of the matrix
 */
class DenseSolver : public BaseMatrixSolver {
private:
    Eigen::MatrixXd matrix;


public:
    // CONSTRUCTORS
    /**
     *  @param matrix                               the full dense representation of the matrix
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
     */
    DenseSolver(const Eigen::MatrixXd& matrix, size_t number_of_requested_eigenpairs = 1);

    /**
     *  Constructor that sets a zero matrix
     *
     *  @param dim                                  the dimension of the matrix
     *  @param number_of_requested_eigenpairs       the number of eigenpairs the eigensolver should find
     */
    DenseSolver(size_t dim, size_t number_of_requested_eigenpairs = 1);

    /**
     *  @param matrix                   the full dense representation of the matrix
     *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
     */
    DenseSolver(const Eigen::MatrixXd& matrix, const DenseSolverOptions& dense_solver_options);

    /**
     *  Constructor that sets a zero matrix
     *
     *  @param dim                      the dimension of the matrix
     *  @param dense_solver_options     the options to be used for the dense eigenproblem algorithm
     */
    DenseSolver(size_t dim, const DenseSolverOptions& dense_solver_options);


    // DESTRUCTOR
    ~DenseSolver() override = default;


    // GETTERS
    const Eigen::MatrixXd& get_matrix() const { return this->matrix; };


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



}  // namespace GQCP



#endif  // GQCP_DENSESOLVER_HPP
