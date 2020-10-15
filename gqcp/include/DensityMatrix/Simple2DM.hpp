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


#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"


namespace GQCP {


/*
 *  MARK: Simple2DM implementation
 */

/**
 *  A two-electron density matrix that is described by a single tensor.
 * 
 *  This class is used as a base class for `Orbital2DM` and `G2DM`, since they are both expressed using a single tensor, as opposed to `SpinResolved2DM`, which uses separate alpha- and beta- tensor. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 *  @tparam _DerivedDM              The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedDM>
class Simple2DM:
    public SquareRankFourTensor<_Scalar> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedDM = _DerivedDM;

    // The type of 'this'.
    using Self = Simple2DM<Scalar, DerivedDM>;

    // The type of the one-electron density matrix that is naturally related to the derived 2-DM.
    using OneDM_Placeholder = typename DensityMatrixTraits<DerivedDM>::OneDM_Placeholder;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit base constructors.
    using SquareRankFourTensor<Scalar>::SquareRankFourTensor;


    /*
     *  MARK: General information
     */
    size_t numberOfOrbitals() const { return this->dimension(); }


    /*
     *  MARK: Contractions
     */

    /**
     *  @return a partial contraction of the 2-DM, where D(p,q) = d(p,q,r,r)
     */
    OneDM_Placeholder reduce() const {

        // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

        const auto K = this->numberOfOrbitals();

        OneDM_Placeholder D = OneDM_Placeholder::Zero(K);
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {

                for (size_t r = 0; r < K; r++) {
                    D(p, q) += this->operator()(p, q, r, r);
                }
            }
        }

        return D;
    }
};


}  // namespace GQCP
