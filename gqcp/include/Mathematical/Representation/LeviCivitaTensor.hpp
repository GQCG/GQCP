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
template <typename _Scalar>
class LeviCivitaTensor {
private:
    // The values of the Levi-Civita tensor.
    Tensor<_Scalar, 3> epsilon;

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The default constructor.
     */
    LeviCivitaTensor();


    /*
     *  MARK: Element access
     */

    /**
     *  @param i            The first accessor index.
     *  @param j            The second accessor index.
     *  @param k            The third accessor index.
     * 
     *  @return The element epsilon(i,j,k) of the Levi-Civita tensor.
     */
    _Scalar operator()(const size_t i, const size_t j, const size_t k) const { return this->epsilon(i, j, k); }

    /**
     *  Find the index k such that the Levi-Civita tensor epsilon(i,j,k) (or any permutations thereof) does not vanish.
     * 
     *  @param i            An index of the Levi-Civita tensor.
     *  @param j            Another index of the Levi-Civita tensor.
     * 
     *  @return The index k such that epsilon(i,j,k) does not vanish.
     */
    size_t nonZeroIndex(const size_t i, const size_t j) const;
};


}  // namespace GQCP
