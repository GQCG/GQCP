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


#include "Basis/Transformations/UOrbitalRotationGeneratorsComponent.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "QuantumChemical/SpinResolved.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"
#include "Utilities/complex.hpp"


namespace GQCP {


/*
 *  MARK: UOrbitalRotationGenerators implementation
 */

/**
 *  A type that encapsulates orbital rotation generators for the alpha- and beta-parts separately.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class UOrbitalRotationGenerators:
    public SpinResolvedBase<UOrbitalRotationGeneratorsComponent<_Scalar>, UOrbitalRotationGenerators<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of 'this'.
    using Self = UOrbitalRotationGenerators<Scalar>;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<UOrbitalRotationGeneratorsComponent<Scalar>, Self>::Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<UOrbitalRotationGeneratorsComponent<Scalar>, UOrbitalRotationGenerators<Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Access
     */

    /**
     *  @return The antisymmetric orbital rotation generator spin resolved matrix kappa.
     */
    const SpinResolved<SquareMatrix<Scalar>> asMatrix() const {

        const auto& alpha_component = this->alpha();
        const auto& beta_component = this->beta();

        const auto kappa_matrix_alpha = SquareMatrix<Scalar>::FromStrictTriangle(alpha_component.asVector());  // Lower triangle only.
        const auto kappa_matrix_beta = SquareMatrix<Scalar>::FromStrictTriangle(beta_component.asVector());    // Lower triangle only.

        const SquareMatrix<Scalar> kappa_matrix_alpha_transpose = kappa_matrix_alpha.transpose();
        const SquareMatrix<Scalar> kappa_matrix_beta_transpose = kappa_matrix_beta.transpose();

        // Add the antisymmetric component and return the matrix representation.
        return SpinResolved<SquareMatrix<Scalar>> {(kappa_matrix_alpha - kappa_matrix_alpha_transpose), (kappa_matrix_beta - kappa_matrix_beta_transpose)};
    }

    /**
     *  @return The orbital rotation generators as the strict upper/lower triangle of the spin resolved kappa matrix.
     */
    const SpinResolved<VectorX<Scalar>>& asVector() const {

        return SpinResolved<VectorX<Scalar>> {this->alpha().asVector(), this->beta().asVector()};
    }


    /**
     *  @return The number of alpha and beta spatial orbitals that can be rotated using these orbital rotation generators.
     */
    SpinResolved<size_t> numberOfSpatialOrbitals() const {

        const auto& alpha_component = this->alpha();
        const auto& beta_component = this->beta();

        return SpinResolved<size_t> {alpha_component.numberOfSpatialOrbitals(), beta_component.numberOfSpatialOrbitals()};
    }
};


}  // namespace GQCP
