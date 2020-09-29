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


#include "Basis/Transformations/SimpleTransformationMatrix.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 *  A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class UTransformationMatrix:
    public SpinResolved<TransformationMatrix<_Scalar>, UTransformationMatrix<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit SpinResolved's constructors.
    using SpinResolved<TransformationMatrix<Scalar>, UTransformationMatrix<Scalar>>::SpinResolved;


    /*
     *  MARK: Named constructors
     */

    /**
     *  @param T                the equal transformation matrix for both the alpha and the beta spin-orbitals
     * 
     *  @return a spin-resolved transformation matrix where the transformation matrices for the alpha and beta spin-orbitals are equal
     */
    static UTransformationMatrix<Scalar> FromRestricted(const TransformationMatrix<Scalar>& T) { return UTransformationMatrix<Scalar>::FromEqual(T); }


    /*
     *  MARK: General information
     */

    /**
     *  @param sigma            Alpha or beta.
     * 
     *  @return The number of orbitals (spin-orbitals) that the transformation matrix for the sigma spin-orbitals is related to
     */
    size_t numberOfOrbitals(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->alpha().numberOfOrbitals();
            break;
        }

        case Spin::beta: {
            return this->beta().numberOfOrbitals();
            break;
        }
        }
    }
};


}  // namespace GQCP
