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


#include "Basis/Transformations/UTransformation.hpp"
#include "Domain/UMullikenDomainComponent.hpp"
#include "QuantumChemical/SpinResolved.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  An unrestricted Mulliken-based domain in an AO basis.
 *
 *  @param _Scalar          The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
 */
template <typename _Scalar>
class UMullikenDomain:
    public SpinResolvedBase<UMullikenDomainComponent<_Scalar>, UMullikenDomain<_Scalar>> {
public:
    // The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
    using Scalar = _Scalar;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<UMullikenDomainComponent<Scalar>, UMullikenDomain<Scalar>>::Of;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<UMullikenDomainComponent<Scalar>, UMullikenDomain<Scalar>>::SpinResolvedBase;


    /**
     *  MARK: Domain and projecting
     */

    /**
     *  @return The partition matrices 'P_A' (alpha and beta) related to this unrestricted Mulliken domain.
     */
    SpinResolved<SquareMatrix<Scalar>> partitionMatrix(const UTransformation<Scalar>& C) const {

        return SpinResolved<SquareMatrix<Scalar>> {this->alpha().partitionMatrix(C.alpha()), this->beta().partitionMatrix(C.beta())};
    }


    /**
     *  @return The Mulliken projection matrix (as an unrestricted transformation) defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix, for both spin components.
     */
    UTransformation<Scalar> projectionMatrix(const UTransformation<Scalar>& C) const {

        return UTransformation<Scalar> {this->alpha().projectionMatrix(C.alpha()), this->beta().projectionMatrix(C.beta())};
    }
};


}  // namespace GQCP
