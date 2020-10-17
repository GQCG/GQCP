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


#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/*
 *  MARK: SimpleTransformationMatrix implementation
 */

/**
 *  A basis transformation that is represented by a single transformation matrix.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 *  
 *  This class is used as a base class for `RTransformationMatrix` and `GTransformationMatrix`, since they are both expressed using a single matrix, as opposed to `UTransformationMatrix`, which uses separate transformation coefficients for alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                                 The scalar type used for a transformation coefficient: real or complex.
 *  @tparam _DerivedTransformationMatrix            The type of the transformation matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedTransformationMatrix>
class SimpleTransformationMatrix:
    public SquareMatrix<_Scalar>,
    public BasisTransformable<_DerivedTransformationMatrix, _DerivedTransformationMatrix>,
    public JacobiRotatable<_DerivedTransformationMatrix> {

public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedTransformationMatrix = _DerivedTransformationMatrix;

    // The type of 'this'.
    using Self = SimpleTransformationMatrix<_Scalar, _DerivedTransformationMatrix>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit SquareMatrix' constructors.
    using SquareMatrix<Scalar>::SquareMatrix;


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a transformation matrix from Jacobi rotation parameters. Note that we work with the (cos, sin, -sin, cos) definition.
     * 
     *  @param jacobi_rotation_parameters               The parameters that define the Jacobi rotation matrix.
     *  @param dim                                      The dimension of the resulting matrix.
     *
     *  @return The Jacobi rotation matrix that corresponds to Jacobi rotation parameters.
     */
    static DerivedTransformationMatrix FromJacobi(const JacobiRotationParameters& jacobi_rotation_parameters, const size_t dim) {

        // Create an identity transformation matrix and apply a Jacobi rotation.
        DerivedTransformationMatrix J = SquareMatrix<Scalar>::Identity(dim);

        return J.rotated(jacobi_rotation_parameters);
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals (spinors, spin-orbitals or spatial orbitals, depending on the context/derived class) this transformation matrix is related to.
     */
    size_t numberOfOrbitals() const { return this->dimension(); /* the dimension of the square matrix */ }


    /*
     *  MARK: Conforming to BasisTransformable
     */

    /**
     *  Apply the basis transformation and return the resulting one-electron integrals.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed one-electron integrals.
     */
    DerivedTransformationMatrix transformed(const DerivedTransformationMatrix& transformation_matrix) const override {

        return DerivedTransformationMatrix {(*this) * transformation_matrix};
    }


    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedTransformationMatrix, DerivedTransformationMatrix>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedTransformationMatrix, DerivedTransformationMatrix>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_parameters        The Jacobi rotation parameters.
     * 
     *  @return The transformation matrix that that encapsulates the sequential application of this transformation, followed by the Jacobi rotation.
     */
    virtual DerivedTransformationMatrix rotated(const JacobiRotationParameters& jacobi_parameters) const override {

        const auto p = jacobi_parameters.p();
        const auto q = jacobi_parameters.q();

        // Create the matrix representation of a Jacobi rotation using Eigen's APIs.
        // We're applying the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T).
        const auto jacobi_rotation = jacobi_parameters.Eigen();
        auto result = (*this);
        result.applyOnTheRight(p, q, jacobi_rotation);

        return DerivedTransformationMatrix {result};
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedTransformationMatrix>::rotate;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename DerivedTransformationMatrix>
struct BasisTransformableTraits<SimpleTransformationMatrix<Scalar, DerivedTransformationMatrix>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix". A transformation matrix should naturally be transformable with itself.
    using TM = DerivedTransformationMatrix;
};


}  // namespace GQCP
