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

#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"

#include <Eigen/Cholesky>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param alpha                the increment factor used to obtain the next scaler tau
 *  @param beta                 the heuristic increment
 */
IterativeIdentitiesHessianModifier::IterativeIdentitiesHessianModifier(const double alpha, const double beta) :
    alpha {alpha},
    beta {beta} {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param hessian      the current indefinite Hessian
 * 
 *  @return a modified Hessian that is made positive (for minimizers) or negative (for maximizers) definite
 */
SquareMatrix<double> IterativeIdentitiesHessianModifier::operator()(const SquareMatrix<double>& hessian) {

    // Initialize tau_0, i.e. the initial value for tau
    const double min_diagonal_el = hessian.diagonal().minCoeff();
    if (min_diagonal_el > 0) {
        this->tau = 0.0;
    } else {
        this->tau = -min_diagonal_el + this->beta;
    }


    // Update the modified hessian with multiples of the identity matrix until it is positive/negative definite
    const size_t dim = hessian.dimension();
    SquareMatrix<double> modified_hessian = hessian;
    Eigen::LLT<Eigen::MatrixXd> llt_factorizer {};  // use Cholesky decomposition to check for positive/negative definiteness

    while ((llt_factorizer.compute(modified_hessian), llt_factorizer.info() != Eigen::Success)) {  // comma operator gets value of last expression

        // Modify the Hessian with a multiple of the identity matrix
        modified_hessian += this->tau * SquareMatrix<double>::Identity(dim, dim);

        // Update the scaler
        this->tau = std::max(this->alpha * this->tau, this->beta);
    }

    return modified_hessian;
}


}  // namespace GQCP
