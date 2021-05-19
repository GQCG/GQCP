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
 *  @tparam _Scalar         The scalar type used for a orbital rotation generator: real or complex.
 */
template <typename _Scalar>
class UOrbitalRotationGenerators:
    public SpinResolvedBase<UOrbitalRotationGeneratorsComponent<_Scalar>, UOrbitalRotationGenerators<_Scalar>> {
public:
    // The scalar type used for a orbital rotation generator: real or complex.
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
     *  @return The anti-Hermitian orbital rotation generator spin resolved matrix kappa.
     */
    const SpinResolved<SquareMatrix<Scalar>> asMatrix() const {

        return SpinResolved<SquareMatrix<Scalar>> {this->alpha().asMatrix(), this->beta().asMatrix()};
    }

    /**
     *  @return The orbital rotation generators as the strict lower triangle of the spin resolved kappa matrix.
     */
    const SpinResolved<VectorX<Scalar>>& asVector() const {

        return SpinResolved<VectorX<Scalar>> {this->alpha().asVector(), this->beta().asVector()};
    }


    /**
     *  @return The number of alpha and beta spatial orbitals that can be rotated using these orbital rotation generators.
     */
    SpinResolved<size_t> numberOfOrbitals() const {

        return SpinResolved<size_t> {this->alpha().numberOfOrbitals(), this->beta().numberOfOrbitals()};
    }
};


}  // namespace GQCP
