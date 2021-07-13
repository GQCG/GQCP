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
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "QuantumChemical/SpinResolved.hpp"

#include <vector>


namespace GQCP {


/**
 *  A generalized Mulliken-based partitioning of an AO basis.
 * 
 *  @param _Scalar          The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
 */
template <typename _Scalar>
class GMullikenPartitioning {
public:
    // The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
    using Scalar = _Scalar;


private:
    // A set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis. They are equal for the alpha- and beta-underlying scalar bases.
    std::vector<size_t> m_indices;

    // The transformation that relates the atomic spinors to the set of current generalized spinors.
    GTransformation<Scalar> C;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a generalized Mulliken partitioning from a set of included AO indices.
     * 
     *  @param indices          A set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis. They are equal for the alpha- and beta-underlying scalar bases.
     *  @param C                The transformation that relates the atomic spinors to the set of current generalized spinors.
     */
    GMullikenPartitioning(const std::vector<size_t>& indices, const GTransformation<Scalar>& C) :
        m_indices {indices},
        C {C} {}


    /**
     *  MARK: General information
     */

    /**
     *  @return A set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis. They are equal for the alpha- and beta-underlying scalar bases.
     */
    const std::vector<size_t>& indices() const { return this->m_indices; }


    /**
     *  @return The number of orbitals that this Mulliken partitioning is related to.
     */
    size_t numberOfOrbitals() const { return this->C.numberOfOrbitals(); }


    /**
     *  MARK: Partitioning and projecting
     */

    /**
     *  @return The partition matrix 'P_A' related to this Mulliken partitioning.
     */
    SquareMatrix<Scalar> partitionMatrix() const {

        const auto M = this->numberOfOrbitals();
        const auto K = M / 2;

        // Set up the top-left (alpha) and bottom-right (beta) blocks of the total partition matrix.
        SquareMatrix<Scalar> P = SquareMatrix<Scalar>::Zero(M);

        auto P_component = SquareMatrix<Scalar>::PartitionMatrix(this->indices(), K);
        P.topLeftCorner(K, K) = P_component;
        P.bottomRightCorner(K, K) = P_component;

        return P;
    }


    /**
     *  @return The Mulliken projection, defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix.
     * 
     *  @note We are aware that this formula is duplicate code (see `RMullikenPartitioning`), but it isn't worth (yet) to refactor this common functionality into a base class.
     */
    GTransformation<Scalar> projectionMatrix() const { return GTransformation<Scalar> {this->C.inverse().matrix() * this->partitionMatrix() * this->C.matrix()}; }
};


}  // namespace GQCP
