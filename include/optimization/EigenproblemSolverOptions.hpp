#include <utility>

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
#ifndef NUMOPT_EIGENPROBLEMSOLVEROPTIONS_HPP
#define NUMOPT_EIGENPROBLEMSOLVEROPTIONS_HPP


#include <cstddef>
#include <Eigen/Dense>


namespace numopt {
namespace eigenproblem {


/**
 *  An enum class for the implemented eigenproblem solver types
 */
enum class SolverType {
    DENSE,
    SPARSE,
    DAVIDSON
};



/**
 *  A base struct to specify eigenproblem solver options, whose derived structs can be used with the eigenproblem solvers
 */
struct BaseSolverOptions {
public:
    // MEMBERS
    size_t number_of_requested_eigenpairs = 1;


    // PURE VIRTUAL METHODS
    virtual numopt::eigenproblem::SolverType get_solver_type() = 0;
};



/**
 *  A struct to specify dense eigenproblem solver options
 */
struct DenseSolverOptions : public BaseSolverOptions {
public:
    // OVERRIDDEN METHODS
    numopt::eigenproblem::SolverType get_solver_type () override { return numopt::eigenproblem::SolverType::DENSE; };
};



/**
 *  A struct to specify sparse eigenproblem solver options
 */
struct SparseSolverOptions : public BaseSolverOptions {
public:
    // OVERRIDDEN METHODS
    numopt::eigenproblem::SolverType get_solver_type () override { return numopt::eigenproblem::SolverType::SPARSE; };
};



/**
 *  A struct to specify Davidson eigenproblem solver options
 */
struct DavidsonSolverOptions : public BaseSolverOptions {
public:
    // MEMBERS
    double residue_tolerance = 1.0e-08;  // the tolerance on the norm of the residual vector
    double correction_threshold = 1.0e-12;  // the threshold used in solving the (approximated) residue correction equation

    size_t maximum_subspace_dimension = 15;
    size_t collapsed_subspace_dimension = 2;

    Eigen::MatrixXd X_0;  // Eigen::MatrixXd of initial guesses, or Eigen::VectorXd of initial guess


    // CONSTRUCTORS
    /**
     *  @param initial_guess        the initial guess(es) for the Davidson algorithm, specified as column(s) of the given vector/matrix
     */
    explicit DavidsonSolverOptions(const Eigen::MatrixXd& initial_guess) :
        X_0 (initial_guess)
    {}


    // OVERRIDDEN METHODS
    numopt::eigenproblem::SolverType get_solver_type () override { return numopt::eigenproblem::SolverType::DAVIDSON; };
};


}  // namespace eigenproblem
}  // namespace numopt



#endif  // NUMOPT_EIGENPROBLEMSOLVEROPTIONS_HPP
