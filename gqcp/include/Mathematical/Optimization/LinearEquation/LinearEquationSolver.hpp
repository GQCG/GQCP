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
#pragma once


#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Optimization/LinearEquation/ColPivHouseholderQRSolution.hpp"
#include "Mathematical/Optimization/LinearEquation/HouseholderQRSolution.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"


namespace GQCP {


/**
 *  A factory class that can construct solvers for linear systems of equations in an easy way.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the variables of the system of equations
 */
template <typename _Scalar>
class LinearEquationSolver {
public:
    using Scalar = _Scalar;


public:
    /*
     * STATIC PUBLIC METHODS
     */

    /**
     *  @return a linear equations solver that uses the Householder QR algorithm
     */
    static Algorithm<LinearEquationEnvironment<Scalar>> HouseholderQR() {

        // Our Householder QR decomposition is just a wrapper around Eigen's.
        StepCollection<LinearEquationEnvironment<Scalar>> householder_steps {};
        householder_steps.add(HouseholderQRSolution<Scalar>());

        return Algorithm<LinearEquationEnvironment<Scalar>>(householder_steps);
    }


    /**
     *  @return a linear equations solver that uses the Householder QR (with column-pivoting) algorithm
     */
    static Algorithm<LinearEquationEnvironment<Scalar>> ColPivHouseholderQR() {

        // Our Householder QR decomposition is just a wrapper around Eigen's.
        StepCollection<LinearEquationEnvironment<Scalar>> householder_steps {};
        householder_steps.add(ColPivHouseholderQRSolution<Scalar>());

        return Algorithm<LinearEquationEnvironment<Scalar>>(householder_steps);
    }
};


}  // namespace GQCP
