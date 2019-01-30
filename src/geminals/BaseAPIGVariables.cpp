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
#include "geminals/BaseAPIGVariables.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param x        the variables in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
BaseAPIGVariables::BaseAPIGVariables(const Eigen::VectorXd& x, size_t N_P, size_t K) :
    N_P (N_P),
    K (K),
    x (x)
{
    // The base constructor just sets the members, derived constructors should perform subsequent checks on the arguments
}


/**
 *  Default constructor setting everything to zero
 */
BaseAPIGVariables::BaseAPIGVariables() :
    BaseAPIGVariables(Eigen::VectorXd::Zero(0), 0, 0)
{}



/*
 *  OPERATORS
 */

/**
 *  @param i        the major (geminal, subscript, non-contiguous) index
 *  @param p        the minor (orbital, superscript, contiguous) index
 *
 *  @return the variable X_i^p
 */
double BaseAPIGVariables::operator()(size_t i, size_t p) const {
    size_t mu = this->vectorIndex(i, p);
    return this->operator()(mu);
}


}  // namespace GQCP
