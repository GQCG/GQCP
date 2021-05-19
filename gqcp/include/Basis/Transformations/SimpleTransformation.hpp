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
#include "Basis/Transformations/OrbitalRotationGeneratorTraits.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"

#include <unsupported/Eigen/MatrixFunctions>


namespace GQCP {


/*
 *  MARK: SimpleTransformation implementation
 */

/**
 *  A basis transformation that can be represented by a single transformation matrix.
 * 
 *  In general, we adopt the convention outlined in (https://gqcg-res.github.io/knowdes/spinor-transformations.html), where the new orbitals' coefficients can be found in the respective **column** of the related transformation matrix.
 *  
 *  This class is used as a base class for `RTransformation` and `GTransformation`, since they are both expressed using a single matrix, as opposed to `UTransformation`, which uses separate transformation coefficients for alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                                 The scalar type used for a transformation coefficient: real or complex.
 *  @tparam _DerivedTransformation                  The type of the transformation matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedTransformation>
class SimpleTransformation:
    public BasisTransformable<_DerivedTransformation>,
    public JacobiRotatable<_DerivedTransformation> {

public:
    // The scalar type used for a transformation coefficient: real or complex.
    using Scalar = _Scalar;

    // The type of the transformation that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedTransformation = _DerivedTransformation;

    // The type of 'this'.
    using Self = SimpleTransformation<Scalar, DerivedTransformation>;

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;

    // The type of the orbital rotation generators that is naturally associated with the derived transformation.
    using OrbitalRotationGeneratorType = typename OrbitalRotationGeneratorTraits<DerivedTransformation>::OrbitalRotationGenerators;


protected:
    // The transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.
    SquareMatrix<Scalar> T;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a `SimpleTransformation` from the transformation matrix that it encapsulates.
     * 
     *  @param T                The transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.
     */
    SimpleTransformation(const SquareMatrix<Scalar>& T) :
        T {T} {}


    /**
     *  Construct a `SimpleTransformation` from the transformation matrix (associated with a set of orbital rotation generators) that it encapsulates.
     * 
     *  @param orbital_rotation_generators                The orbital rotation generators from which a transformation so-called `kappa matrix` is constructed.
     */
    SimpleTransformation(const OrbitalRotationGeneratorType& orbital_rotation_generators) :
        T {(-orbital_rotation_generators.asMatrix()).exp()} {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a general transformation from Jacobi rotation. Note that we work with the (cos, sin, -sin, cos) definition.
     * 
     *  @param jacobi_rotation                      The Jacobi rotation.
     *  @param dim                                  The dimension of the resulting matrix.
     *
     *  @return The general transformation that corresponds to the given Jacobi rotation.
     */
    static DerivedTransformation FromJacobi(const JacobiRotation& jacobi_rotation, const size_t dim) {

        // Create an identity transformation matrix and apply a Jacobi rotation.
        DerivedTransformation J = SquareMatrix<Scalar>::Identity(dim);

        return J.rotated(jacobi_rotation);
    }


    /**
     *  Create an identity transformation between two orbital bases.
     * 
     *  @param dim              The dimension of the transformation matrix.
     */
    static DerivedTransformation Identity(const size_t dim) { return DerivedTransformation {SquareMatrix<Scalar>::Identity(dim)}; }

    /**
     *  Create a random transformation.
     * 
     *  @param dim          The dimension of the transformation matrix.
     */
    static DerivedTransformation Random(const size_t dim) { return DerivedTransformation {SquareMatrix<Scalar>::Random(dim)}; }

    /**
     *  Create a random unitary transformation.
     * 
     *  @param dim          The dimension of the transformation matrix.
     */
    static DerivedTransformation RandomUnitary(const size_t dim) { return DerivedTransformation {SquareMatrix<Scalar>::RandomUnitary(dim)}; }

    /**
     *  Create a zero transformation.
     * 
     *  @param dim              The dimension of the transformation matrix.
     */
    static DerivedTransformation Zero(const size_t dim) { return DerivedTransformation {SquareMatrix<Scalar>::Zero(dim)}; }


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of orbitals (spinors, spin-orbitals or spatial orbitals, depending on the context/derived class) this transformation is related to.
     */
    size_t numberOfOrbitals() const { return this->matrix().dimension(); }

    /**
     *  @return The dimension of this basis transformation.
     */
    size_t dimension() const { return this->numberOfOrbitals(); }


    /*
     *  MARK: Transformation matrix
     */

    /**
     *  @return The transformation matrix that collects the expansion coefficients of the new basis (vectors) in the old basis as columns.
     */
    const SquareMatrix<Scalar>& matrix() const { return this->T; }


    /*
     *  MARK: Linear algebra
     */

    /**
     *  @return The adjoint transformation of this one.
     */
    DerivedTransformation adjoint() const { return DerivedTransformation {this->matrix().adjoint()}; }

    /**
     *  @return The inverse transformation of this one.
     */
    DerivedTransformation inverse() const { return DerivedTransformation {this->matrix().inverse()}; }

    /**
     *  Check if this transformation is unitary.
     * 
     *  @param threshold                The threshold used to check for unitarity.
     * 
     *  @return If this transformation is unitary, within the given threshold.
     */
    bool isUnitary(const double threshold = 1.0e-12) const { return this->matrix().isUnitary(threshold); }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result, which corresponds to the concatenation of two basis transformations.
     * 
     *  @param T        The basis transformation.
     * 
     *  @return The transformation that encapsulates the sequential application of this transformation, followed by the given transformation.
     */
    DerivedTransformation transformed(const DerivedTransformation& T) const override { return DerivedTransformation {this->matrix() * T.matrix()}; }

    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedTransformation>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedTransformation>::rotated;


    /*
     *  MARK: Conforming to JacobiRotatable
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The transformation that encapsulates the sequential application of this transformation, followed by the Jacobi rotation.
     */
    DerivedTransformation rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto p = jacobi_rotation.p();
        const auto q = jacobi_rotation.q();

        // Create the matrix representation of a Jacobi rotation using Eigen's APIs.
        // We're applying the Jacobi rotation as J = I * jacobi_rotation (cfr. B' = B T).
        const auto eigen_jacobi_rotation = jacobi_rotation.Eigen();
        auto result = this->matrix();
        result.applyOnTheRight(p, q, eigen_jacobi_rotation);

        return DerivedTransformation {result};
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedTransformation>::rotate;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename DerivedTransformation>
struct BasisTransformableTraits<SimpleTransformation<Scalar, DerivedTransformation>> {

    // The type of the transformation for which the basis transformation should be defined. A transformation matrix should naturally be transformable with itself.
    using Transformation = DerivedTransformation;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 * 
 *  @tparam T       The type that should conform to `JacobiRotatable`.
 */
template <typename Scalar, typename DerivedTransformation>
struct JacobiRotatableTraits<SimpleTransformation<Scalar, DerivedTransformation>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
