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


#include "Basis/Transformations/RTransformationMatrix.hpp"
#include "Basis/Transformations/UTransformationMatrixComponent.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/*
 *  MARK: UTransformationMatrix implementation
 */

/**
 *  A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases.
 * 
 *  @tparam _Scalar         The scalar type used for a transformation coefficient: real or complex.
 */
template <typename _Scalar>
class UTransformationMatrix:
    public SpinResolvedBase<UTransformationMatrixComponent<_Scalar>, UTransformationMatrix<_Scalar>>,
    public BasisTransformable<UTransformationMatrix<_Scalar>> {
public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation matrix for which the basis transformation should be defined. A transformation matrix should naturally be transformable with itself.
    using TM = UTransformationMatrix<Scalar>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<UTransformationMatrixComponent<Scalar>, UTransformationMatrix<Scalar>>::SpinResolvedBase;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create an UTransformationMatrix from an RTransformationMatrix, leading to transformation matrices for both spin components that are equal.
     * 
     *  @param T                The equal transformation matrix for both the alpha and the beta spin-orbitals.
     * 
     *  @return An UTransformationMatrix.
     */
    static UTransformationMatrix<Scalar> FromRestricted(const RTransformationMatrix<Scalar>& T) { return UTransformationMatrix<Scalar>::FromEqual(T); }


    /**
     *  Create an identity UTransformationMatrix.
     * 
     *  @param dim_alpha            The number of alpha spin-orbitals.
     *  @param dim_beta             The number of beta spin-orbitals.
     * 
     *  @return An identity UTransformationMatrix.
     */
    static UTransformationMatrix<Scalar> Identity(const size_t dim_alpha, const size_t dim_beta) {

        const UTransformationMatrixComponent<Scalar> T_alpha = UTransformationMatrixComponent<Scalar>::Identity(dim_alpha);
        const UTransformationMatrixComponent<Scalar> T_beta = UTransformationMatrixComponent<Scalar>::Identity(dim_beta);

        return UTransformationMatrix<Scalar> {T_alpha, T_beta};
    }


    /**
     *  Create an identity UTransformationMatrix.
     * 
     *  @param dim              The dimension of the alpha and beta spin-orbitals.
     * 
     *  @return An identity UTransformationMatrix.
     */
    static UTransformationMatrix<Scalar> Identity(const size_t dim) { return UTransformationMatrix<Scalar>::Identity(dim, dim); }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    UTransformationMatrix<Scalar> transformed(const TM& transformation_matrix) const override {

        // Transform the components of 'this' with the components of the transformation matrix.
        auto result = *this;
        result.alpha().transform(transformation_matrix.alpha());
        result.beta().transform(transformation_matrix.beta());

        return result;
    }


    /*
     *  MARK: General information
     */

    /**
     *  @param threshold        The threshold used for checking unitarity.
     * 
     *  @return If this transformation matrix is considered to be unitary, within the given treshold.
     */
    bool isUnitary(const double threshold) const { return this->alpha().isUnitary(threshold) && this->beta().isUnitary(threshold); }
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar>
struct BasisTransformableTraits<UTransformationMatrix<Scalar>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix". A transformation matrix should naturally be transformable with itself.
    using TM = UTransformationMatrix<Scalar>;
};


}  // namespace GQCP
