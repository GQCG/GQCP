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


#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Domain/DiscreteDomain.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 * A Mulliken partitioned domain as a collection of atomic orbitals.
 * The sites {`i`} that are present in the domain are represented by a set bit at the corresponding indices `i`.
 */
template <typename _Scalar>
class UMullikenDomainComponent:
    public DiscreteDomain {
public:
    using DiscreteDomain::DiscreteDomain;

    // The expansion scalar type used for a Mulliken domain: real or complex.
    using Scalar = _Scalar;

    /**
     *  MARK: domain and projecting
     */

    /**
     *  @return The partition matrix 'P_A' related to this Mulliken domain.
     */
    SquareMatrix<Scalar> partitionMatrix(const UTransformationComponent<Scalar>& C) const { return SquareMatrix<Scalar>::PartitionMatrix(this->domainIndices(), C.numberOfOrbitals()); }


    /**
     * @param C         The transformation that relates the atomic spinors to the set of current spin orbitals.
     *
     *  @return The Mulliken projection, defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix.
     *
     *  @note We are aware that this formula is duplicate code (see `GMullikenDomain`), but it isn't worth (yet) to refactor this common functionality into a base class.
     */
    UTransformationComponent<Scalar> projectionMatrix(const UTransformationComponent<Scalar>& C) const { return UTransformationComponent<Scalar> {C.inverse().matrix() * this->partitionMatrix(C) * C.matrix()}; }
};


}  // namespace GQCP
