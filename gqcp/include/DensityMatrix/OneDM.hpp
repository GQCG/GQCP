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


#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A type that represents a one-electron density matrix.
 *
 *  @tparam _Scalar     the scalar type
 */
template <typename _Scalar>
class OneDM: public SquareMatrix<_Scalar> {
public:
    using Scalar = _Scalar;
    using Self = OneDM<Scalar>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit base constructors.
    using SquareMatrix<Scalar>::SquareMatrix;


    /*
     *  MARK: General information
     */
    size_t numberOfOrbitals() const { return this->dimension(); }


    /*
     *  MARK: Transformations
     */

    /**
     *  @param T          transformation matrix for the spin unresolved 1-DM
     * 
     *  @return the transformed density matrix.
     */
    OneDM<Scalar> transformed(const TransformationMatrix<double>& T) const {
        return Self(T.adjoint() * (*this) * T);
    }
};


}  // namespace GQCP
