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


#include "Mathematical/Optimization/Minimization/BaseHessianModifier.hpp"


namespace GQCP {


/**
 *  A Hessian modifier functor that uses multiples of the identity matrix to try to obtain a positive/negative definite Hessian
 * 
 *  This implements Algorithm 3.3 in Nocedal and Wright
 */
class IterativeIdentitiesHessianModifier: public BaseHessianModifier {
private:
    double alpha;  // the scaling factor used to obtain the next scaler tau
    double beta;   // the heuristic increment
    double tau;    // the scaling scalar (the 'scaler')


public:
    // CONSTRUCTORS

    /**
     *  @param alpha                the increment factor used to obtain the next scaler tau
     *  @param beta                 the heuristic increment
     */
    IterativeIdentitiesHessianModifier(const double alpha = 2.0, const double beta = 1.0e-03);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param hessian          the current indefinite Hessian
     * 
     *  @return a modified Hessian that is made positive (for minimizers) or negative (for maximizers) definite
     */
    SquareMatrix<double> operator()(const SquareMatrix<double>& hessian) override;
};


}  // namespace GQCP
