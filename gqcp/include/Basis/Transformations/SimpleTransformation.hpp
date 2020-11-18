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
 *  MARK: SimpleTransformation implementation
 */

/**
 *  A basis transformation that is represented by a single transformation matrix.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 *  
 *  This class is used as a base class for `RTransformation` and `GTransformation`, since they are both expressed using a single matrix, as opposed to `UTransformation`, which uses separate transformation coefficients for alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                                 The scalar type used for a transformation coefficient: real or complex.
 *  @tparam _DerivedTransformationMatrix            The type of the transformation matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedTransformationMatrix>
class SimpleTransformation:
    public SquareMatrix<_Scalar>,
    public BasisTransformable<_DerivedTransformationMatrix>,
    public JacobiRotatable<_DerivedTransformationMatrix> {

public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedTransformationMatrix = _DerivedTransformationMatrix;

    // The type of 'this'.
    using Self = SimpleTransformation<_Scalar, _DerivedTransformationMatrix>;

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;


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
     *  Create a transformation matrix from Jacobi rotation. Note that we work with the (cos, sin, -sin, cos) definition.
     * 
     *  @param jacobi_rotation                      The Jacobi rotation.
     *  @param dim                                  The dimension of the resulting matrix.
     *
     *  @return The Jacobi rotation matrix that corresponds to Jacobi rotation.
     */
    static DerivedTransformationMatrix FromJacobi(const JacobiRotation& jacobi_rotation, const size_t dim) {

        // Create an identity transformation matrix and apply a Jacobi rotation.
        DerivedTransformationMatrix J = SquareMatrix<Scalar>::Identity(dim);

        return J.rotated(jacobi_rotation);
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
    using BasisTransformable<DerivedTransformationMatrix>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedTransformationMatrix>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The transformation matrix that that encapsulates the sequential application of this transformation, followed by the Jacobi rotation.
     */
    DerivedTransformationMatrix rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto p = jacobi_rotation.p();
        const auto q = jacobi_rotation.q();

        // Create the matrix representation of a Jacobi rotation using Eigen's APIs.
        // We're applying the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T).
        const auto eigen_jacobi_rotation = jacobi_rotation.Eigen();
        auto result = (*this);
        result.applyOnTheRight(p, q, eigen_jacobi_rotation);

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
struct BasisTransformableTraits<SimpleTransformation<Scalar, DerivedTransformationMatrix>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix". A transformation matrix should naturally be transformable with itself.
    using TM = DerivedTransformationMatrix;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 * 
 *  @tparam T       The type that should conform to `JacobiRotatable`.
 */
template <typename Scalar, typename DerivedTransformationMatrix>
struct JacobiRotatableTraits<SimpleTransformation<Scalar, DerivedTransformationMatrix>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
