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
#include "WaveFunction/WaveFunction.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  @param base_fock_space      the Fock space in which the wave function 'lives'
 *  @param coefficients         the expansion coefficients
 */
WaveFunction::WaveFunction(const BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients) :
    fock_space (BaseFockSpace::CloneToHeap(base_fock_space)),
    coefficients (coefficients)
{}




/*
 *  PUBLIC METHODS
 */
/**
 *  @return the Shannon entropy (or information content) of the wave function
 */
double WaveFunction::calculateShannonEntropy() const {

    // Sum over the Fock space dimension, and only include the term if c_k != 0
    // We might as well replace all coeffients that are 0 by 1, since log(1) = 0 so there is no influence on the final entropy value
    Eigen::ArrayXd coefficients_replaced = this->coefficients.unaryExpr([](double c) { return c < 1.0e-18 ? 1 : c;});  // replace 0 by 1

    Eigen::ArrayXd coefficients_squared = coefficients_replaced.square();
    Eigen::ArrayXd log_coefficients_squared = coefficients_squared.log();  // natural logarithm (ln)

    return - 1 / std::log(2) * (coefficients_squared * log_coefficients_squared).sum();
}


}  // namespace GQCP
