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
 *  A class that represents a normalized basis function and can be evaluated at a point in Euclidean space
 */
class BasisFunction : public LinearCombination<double, CartesianGTO> {
private:
    double N_total;  // the total normalization factor


public:
    using Base = LinearCombination<double, CartesianGTO>;


public:
    // CONSTRUCTORS
    /**
     *  @param lc       a linear combination of CartesianGTOs
     */
    BasisFunction(const Base& lc);


    // PUBLIC METHODS
    /**
     *  @return the total normalization factor of this basis function
     */
    double calculateNormalizationFactor() const;  // TODO: implement
};


}  // namespace GQCP


#endif  /* BasisFunction_hpp */
