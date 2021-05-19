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


#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  MARK: SimpleOrbitalRotationGenerator implementation
 */


/**
 *  A set of orbital rotation generators that can be represented by a single vector.
 *   
 *  This class is used as a base class for `ROrbitalRotationGenerator` and `GOrbitalRotationGenerator`, since they are both expressed using a single vector, as opposed to `UOrbitalRotationGenerator`, which uses separate kappa vectors for alpha- and beta- generators. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                                           The scalar type used for a orbital rotation generator vector: real or complex.
 *  @tparam _DerivedOrbitalRotationGenerator                  The type of the orbital rotation generator that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedOrbitalRotationGenerators>
class SimpleOrbitalRotationGenerators {

public:
    // The scalar type used for a orbital rotation generator vector: real or complex.
    using Scalar = _Scalar;

    // The type of the orbital rotation generator that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedOrbitalRotationGenerators = _DerivedOrbitalRotationGenerators;

private:
    // The number of spatial orbitals that can be rotated using these orbital rotation generators.
    size_t number_of_spatial_orbitals;

    // The strict lower triangle of the kappa matrix.
    VectorX<Scalar> v;


public:
    /**
     *  MARK: Constructors
     */

    /**
     *  Create a `SimpleOrbitalRotationGenerator' from a given kappa vector.
     * 
     *  @param  kappa_vector        The orbital rotation generators represented as a vector that corresponds to the strict upper/lower triangle of the kappa matrix.
     */
    SimpleOrbitalRotationGenerators(const VectorX<Scalar>& kappa_vector) :
        kappa_vector {kappa_vector},
        number_of_spatial_orbitals {strictTriangularRootOf(kappa_vector.size())} {}


    /**
     *  Create a `SimpleOrbitalRotationGenerator' from a given kappa matrix.
     * 
     *  @param  kappa_matrix        The orbital rotation generators represented as a vector that corresponds to the full antisymmetric the kappa matrix.
     */
    SimpleOrbitalRotationGenerators(const SquareMatrix<Scalar>& kappa_matrix) :
        SimpleOrbitalRotationGenerators(kappa_matrix.pairWiseStrictReduced()) {}


    /**
     *  MARK: Named constructors
     */

    /**
     *  Construct orbital rotation generators by adding redundant (i.e. 0) occupied-virtual and virtual-virtual generators to the given occupied-occupied generators.
     * 
     *  @param occ_occ_generators       The occupied-occupied orbital rotation generators.
     *  @param K                    The total number of spatial orbitals.
     * 
     *  @return The 'full' orbital rotation generators from the given occupied-occupied generators.
     */
    static DerivedOrbitalRotationGenerators FromOccOcc(const DerivedOrbitalRotationGenerators& o_o_generators, const size_t K) {
        SquareMatrix<Scalar> kappa_full_matrix = SquareMatrix<Scalar>::Zero(K);

        kappa_full_matrix.topLeftCorner(o_o_generators.numberOfSpatialOrbitals(), o_o_generators.numberOfSpatialOrbitals()) = o_o_generators.asMatrix();
        return DerivedOrbitalRotationGenerators(kappa_full_matrix);
    }


    /**
     *  MARK: Access
     */

    /**
     *  @return The antisymmetric orbital rotation generator matrix kappa.
     */
    const SquareMatrix<Scalar> asMatrix() const {

        const auto kappa_matrix = SquareMatrix<Scalar>::FromStrictTriangle(this->kappa_vector);  // Lower triangle only.
        const SquareMatrix<Scalar> kappa_matrix_transpose = kappa_matrix.transpose();

        // Add the antisymmetric component and return the matrix representation.
        return kappa_matrix - kappa_matrix_transpose;
    }

    /**
     *  @return The orbital rotation generators as the strict upper/lower triangle of the kappa matrix.
     */
    const VectorX<Scalar>& asVector() const { return this->kappa_vector; }


    /**
     *  @return The number of spatial orbitals that can be rotated using these orbital rotation generators.
     */
    size_t numberOfSpatialOrbitals() const { return this->number_of_spatial_orbitals; }
};


}  // namespace GQCP
