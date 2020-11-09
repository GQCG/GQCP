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


#include "Basis/Transformations/UTransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "QuantumChemical/SpinResolved.hpp"

#include <vector>


namespace GQCP {


/**
 *  An unrestricted Mulliken-based partitioning of an AO basis.
 * 
 *  @param _Scalar          The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
 */
template <typename _Scalar>
class UMullikenPartitioning {
public:
    // The scalar type used to represent an element of the Mulliken projection matrix: real or complex.
    using Scalar = _Scalar;


private:
    // A set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis.
    std::vector<size_t> m_indices;

    // The transformation that relates the atomic spin-orbitals to the set of current unrestricted spin-orbitals.
    UTransformationMatrix<Scalar> C;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create an unrestricted Mulliken partitioning from a set of included AO indices.
     * 
     *  @param indices          A set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis. They are equal for the alpha- and beta-underlying scalar bases.
     *  @param C                The transformation that relates the atomic spin-orbitals to the set of current unrestricted spin-orbitals.
     */
    UMullikenPartitioning(const std::vector<size_t>& indices, const UTransformationMatrix<Scalar>& C) :
        m_indices {indices},
        C {C} {}


    /**
     *  MARK: General information
     */

    /**
     *  @return The set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis.
     */
    const std::vector<size_t>& indices() const { return this->m_indices; }


    /**
     *  MARK: Partitioning and projecting
     */

    /**
     *  @return The partition matrices 'P_A' (alpha and beta) related to this Mulliken partitioning.
     */
    SpinResolved<SquareMatrix<Scalar>> partitionMatrix() const {

        const auto P = SquareMatrix<Scalar>::PartitionMatrix(this->indices(), this->indices().size());
        return SpinResolved<SquareMatrix<Scalar>>::FromEqual(P);
    }


    /**
     *  @return The Mulliken projection matrix defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix.
     */
    UTransformationMatrix<Scalar> projectionMatrix() const {

        const UTransformationMatrixComponent<Scalar> P_alpha = this->C.alpha().inverse() * this->partitionMatrix().alpha() * this->C.alpha();

        const UTransformationMatrixComponent<Scalar> P_beta = this->C.beta().inverse() * this->partitionMatrix().beta() * this->C.beta();

        UTransformationMatrix<Scalar> {P_alpha, P_beta};
    }
};


}  // namespace GQCP
