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
 *  @tparam _Scalar         The scalar type used for an orbital rotation generator: real or complex.
 */
template <typename _Scalar>
class UOrbitalRotationGenerators:
    public SpinResolvedBase<UOrbitalRotationGeneratorsComponent<_Scalar>, UOrbitalRotationGenerators<_Scalar>> {
public:
    // The scalar type used for an orbital rotation generator: real or complex.
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
     * MARK Named constructors
     */

    /**
     *  Construct orbital rotation generators by adding redundant (i.e. 0) generators to the given occupation_type - occupation_type generators. This is done for the alpha and beta component separately.
     * 
     *  @param generators                           The orbital rotation generators of the specified ocupation types.
     *  @param row_occupation_type                  The occupation type of the rows of the orbital rotation generator kappa matrix.
     *  @param column_occupation_type               The occupation type of the column of the orbital rotation generator kappa matrix.
     *  @param K                                    The total number of orbitals. In unrestricted and unrestricted these will be spin-orbitals.
     * 
     *  @return The 'full' orbital rotation generators from the given row_occupation_type - column_occupation_type generators, separately constructed for the alpha and beta component.
     */
    static UOrbitalRotationGenerators<Scalar> FromOccupationTypes(const UOrbitalRotationGenerators& generators, const OccupationType row_occupation_type, const OccupationType column_occupation_type, const size_t K) {

        return UOrbitalRotationGenerators<Scalar> {UOrbitalRotationGeneratorsComponent<Scalar>::FromOccupationTypes(generators.alpha(), row_occupation_type, column_occupation_type, K), UOrbitalRotationGeneratorsComponent<Scalar>::FromOccupationTypes(generators.beta(), row_occupation_type, column_occupation_type, K)};
    }


    /*
     *  MARK: Access
     */

    /**
     *  @return These orbital rotation generators as spin-resolved anti-Hermitian matrices.
     */
    const SpinResolved<SquareMatrix<Scalar>> asMatrix() const {

        return SpinResolved<SquareMatrix<Scalar>> {this->alpha().asMatrix(), this->beta().asMatrix()};
    }

    /**
     *  @return The orbital rotation generators as the strict lower triangle of the spin resolved kappa matrix.
     */
    SpinResolved<VectorX<Scalar>> asVector() const {

        return SpinResolved<VectorX<Scalar>> {this->alpha().asVector(), this->beta().asVector()};
    }


    /**
     *  @return The number of alpha and beta spin-orbitals that can be rotated using these orbital rotation generators.
     */
    SpinResolved<size_t> numberOfOrbitals() const {

        return SpinResolved<size_t> {this->alpha().numberOfOrbitals(), this->beta().numberOfOrbitals()};
    }
};


}  // namespace GQCP
