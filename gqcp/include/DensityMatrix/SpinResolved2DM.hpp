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


#include "Basis/Transformations/DoublySpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/DoublySpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UJacobiRotation.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "DensityMatrix/MixedSpinResolved2DMComponent.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "DensityMatrix/PureSpinResolved2DMComponent.hpp"
#include "QuantumChemical/DoublySpinResolvedBase.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  A type that encapsulates the spin parts of the spin-resolved two-electron density matrix.
 *
 *  @tparam _Scalar         The scalar type of one of the density matrix elements: real or complex.
 */
template <typename _Scalar>
class SpinResolved2DM:
    public DoublySpinResolvedBase<PureSpinResolved2DMComponent<_Scalar>, MixedSpinResolved2DMComponent<_Scalar>, SpinResolved2DM<_Scalar>>,
    public DoublySpinResolvedBasisTransformable<SpinResolved2DM<_Scalar>>,
    public DoublySpinResolvedJacobiRotatable<SpinResolved2DM<_Scalar>> {
public:
    // The scalar type of one of the density matrix elements: real or complex.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = SpinResolved2DM<Scalar>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `DoublySpinResolvedBase`'s constructors.
    using DoublySpinResolvedBase<PureSpinResolved2DMComponent<_Scalar>, MixedSpinResolved2DMComponent<_Scalar>, SpinResolved2DM<_Scalar>>::DoublySpinResolvedBase;


    /*
     *  MARK: General information
     */

    /**
     *  @param sigma            Alpha or beta.
     *  @param tau              Alpha or beta.
     * 
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the spin-resolved 2-DM.
     */
    size_t numberOfOrbitals(const Spin sigma, const Spin tau) const {

        if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaAlpha().numberOfOrbitals();
        } else if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaBeta().numberOfOrbitals();
        } else if (sigma == Spin::beta && tau == Spin::alpha) {
            return this->betaAlpha().numberOfOrbitals();
        } else {
            return this->betaBeta().numberOfOrbitals();
        }
    }


    /**
     *  @return The orbital (total, spin-summed) two-electron density matrix.
     */
    Orbital2DM<Scalar> orbitalDensity() const {
        return Orbital2DM<Scalar> {this->alphaAlpha().Eigen() + this->alphaBeta().Eigen() + this->betaAlpha().Eigen() + this->betaBeta().Eigen()};
    }


    /**
     *  MARK: Enabling basis transformations
     */

    // Since `rotate` and `rotated` are both defined in `DoublySpinResolvedBasisTransformable` and `DoublySpinResolvedJacobiRotatable`, we have to explicitly enable these methods here.

    // Allow the `rotate` method from `DoublySpinResolvedBasisTransformable`, since there's also a `rotate` from `DoublySpinResolvedJacobiRotatable`.
    using DoublySpinResolvedBasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `DoublySpinResolvedBasisTransformable`, since there's also a `rotated` from `DoublySpinResolvedJacobiRotatable`.
    using DoublySpinResolvedBasisTransformable<Self>::rotated;

    // Allow the `rotate` method from `DoublySpinResolvedJacobiRotatable`, since there's also a `rotate` from `DoublySpinResolvedBasisTransformable`.
    using DoublySpinResolvedJacobiRotatable<Self>::rotate;

    // Allow the `rotated` method from `DoublySpinResolvedJacobiRotatable`, since there's also a `rotated` from `DoublySpinResolvedBasisTransformable`.
    using DoublySpinResolvedJacobiRotatable<Self>::rotated;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<SpinResolved2DM<Scalar>> {

    // The type of transformation that is naturally related to a `SpinResolved2DM`.
    using Transformation = UTransformation<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar>
struct JacobiRotatableTraits<SpinResolved2DM<Scalar>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
