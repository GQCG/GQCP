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


#include "Mathematical/Representation/SquareRankFourTensor.hpp"


namespace GQCP {


/**
 *  One of the pure (i.e. alpha-alpha or beta-beta) spin components of a spin-resolved 2-DM.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 */
template <typename _Scalar>
class PureSpinResolved2DMComponent:
    public SquareRankFourTensor<_Scalar> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SquareRankFourTensor`'s constructors.
    using SquareRankFourTensor<Scalar>::SquareRankFourTensor;


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals that are related to this 2-DM.
     */
    size_t numberOfOrbitals() const { return this->dimension(); }
};


}  // namespace GQCP
