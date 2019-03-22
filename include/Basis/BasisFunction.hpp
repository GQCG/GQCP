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
#ifndef BasisFunction_hpp
#define BasisFunction_hpp


#include "Basis/CartesianGTO.hpp"
#include "math/LinearCombination.hpp"



namespace GQCP {


/**
 *  A class that represents a normalized basis function: it is a contraction of CartesianGTOs with the same Cartesian exponents (i.e. only s-type, px-type, ... contractions are allowed)
 *
 *
 *  The difference with a Shell is that a Shell represents (possibly) multiple basis functions and a BasisFunction can be evaluated on a point in Euclidian space
 */
class BasisFunction : public LinearCombination<double, CartesianGTO> {
private:
    double N_total;  // the total normalization factor


public:
    using Base = LinearCombination<double, CartesianGTO>;


public:
    // CONSTRUCTORS
    /**
     *  Construct a BasisFunction from its underlying linear combination of CartesianGTOs, i.e. its underlying contraction
     *
     *  @param lc       a linear combination of CartesianGTOs, a contraction of CartesianGTOs
     *
     * Note that:
     *      - the CartesianGTOs should all have the same Cartesian exponents
     *      - in order to behave like a normalized linear combination, the coefficients are all multiplied by the total normalization factor
     */
    BasisFunction(const Base& lc);


    // PUBLIC METHODS
    /**
     *  @return the total normalization factor of this basis function
     */
    double calculateNormalizationFactor() const;
};


}  // namespace GQCP


#endif  /* BasisFunction_hpp */
