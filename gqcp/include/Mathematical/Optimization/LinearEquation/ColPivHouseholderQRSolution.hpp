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

#pragma once


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A step that calculates the solution vector x (in Ax=b) through a Housholder QR (with column-pivoting) decomposition.
 * 
 *  @tparam _Scalar             the scalar type of the elements of the vectors and matrices
 */
template <typename _Scalar>
class ColPivHouseholderQRSolution:
    public Step<LinearEquationEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Read the left- and right-hand sides of the system of equations from the environment, and write the solution back to the environment.";
    }


    /**
     *  Read the left- and right-hand sides of the system of equations from the environment, and write the solution back to the environment.
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(LinearEquationEnvironment<Scalar>& environment) override {

        const auto& A = environment.A;  // the left-hand side
        const auto& b = environment.b;  // the right-hand side

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::ColPivHouseholderQR<MatrixType> linear_solver(A);  // this does the Householder QR (with column-pivoting) decomposition

        environment.x = linear_solver.solve(b);
    }
};


}  // namespace GQCP
