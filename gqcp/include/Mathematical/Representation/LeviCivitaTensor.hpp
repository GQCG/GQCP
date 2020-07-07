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


#include "Mathematical/Representation/Tensor.hpp"


namespace GQCP {


/**
 *  The antisymmetric rank-three Levi-Civita tensor.
 */
class LeviCivitaTensor {
private:
    Tensor<double, 3> epsilon;  // the values of the Levi-Civita tensor

public:
    // CONSTRUCTORS

    /**
     *  The default constructor.
     */
    LeviCivitaTensor();


    // OPERATORS

    /**
     *  @param i            the first accessor index
     *  @param j            the second accessor index
     *  @param k            the third accessor index
     * 
     *  @return the element epsilon(i,j,k) of the Levi-Civita tensor.
     */
    double operator()(const size_t i, const size_t j, const size_t k) const { return this->epsilon(i, j, k); }


    // PUBLIC METHODS

    /**
     *  Find the index k such that the Levi-Civita tensor epsilon(i,j,k) (or any permutations thereof) does not vanish.
     * 
     *  @param i            an index of the Levi-Civita tensor
     *  @param j            another index of the Levi-Civita tensor
     * 
     *  @return the index k such that epsilon(i,j,k) does not vanish
     */
    size_t nonZeroIndex(const size_t i, const size_t j) const;
};


}  // namespace GQCP
