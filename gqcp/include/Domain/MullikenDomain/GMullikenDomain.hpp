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


#include "Basis/Transformations/GTransformation.hpp"
#include "Domain/DiscreteDomain.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 * A Mulliken partitioned domain as a collection of atomic orbitals.
 * The sites {`i'} that are present in the domain are represented by a set bit at the corresponding indices `i'.
 */
template <typename _Scalar>
class GMullikenDomain:
    public DiscreteDomain {
public:
    using DiscreteDomain::DiscreteDomain;

    // The expansion scalar type used for a Mulliken domain: real or complex.
    using Scalar = _Scalar;

    /**
     *  MARK: domain and projecting
     */

    /**
     * @param C     The transformation that relates the atomic spinors to the set of current generalized spinors.
     *
     *  @return The partition matrix 'P_A' related to this Mulliken domain.
     */
    SquareMatrix<Scalar> partitionMatrix(const GTransformation<Scalar>& C) const {

        const auto M = C.numberOfOrbitals();
        const auto K = M / 2;

        // Set up the top-left (alpha) and bottom-right (beta) blocks of the total partition matrix.
        SquareMatrix<Scalar> P = SquareMatrix<Scalar>::Zero(M);

        auto P_component = SquareMatrix<Scalar>::PartitionMatrix(this->domain_indices, K);
        P.topLeftCorner(K, K) = P_component;
        P.bottomRightCorner(K, K) = P_component;

        return P;
    }


    /**
     * @param C         The transformation that relates the atomic spinors to the set of current generalized spinors.
     *
     *  @return The Mulliken projection, defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix.
     *
     *  @note We are aware that this formula is duplicate code (see `RMullikenDomain`), but it isn't worth (yet) to refactor this common functionality into a base class.
     */
    GTransformation<Scalar> projectionMatrix(const GTransformation<Scalar>& C) const { return GTransformation<Scalar> {C.inverse().matrix() * this->partitionMatrix(C) * C.matrix()}; }
};
}  // namespace GQCP
