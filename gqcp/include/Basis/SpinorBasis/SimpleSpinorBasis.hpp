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
#include "Utilities/CRTP.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {

/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename SpinorBasis>
struct SpinorBasisTraits {};


/*
 *  MARK: SimpleSpinorBasis
 */

/**
 *  A class that represents a spinor basis that has no internal structure (hence 'simple') with respect to spin components.
 * 
 *  @tparam _ExpansionScalar            The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
 *  @tparam _FinalSpinorBasis           The spinor basis that ultimately derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _ExpansionScalar, typename _FinalSpinorBasis>
class SimpleSpinorBasis:
    public CRTP<_FinalSpinorBasis>,
    public BasisTransformable<_FinalSpinorBasis>,
    public JacobiRotatable<_FinalSpinorBasis> {

public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The spinor basis that ultimately derives from this class, enabling CRTP and compile-time polymorphism.
    using FinalSpinorBasis = _FinalSpinorBasis;

    // The type of transformation that is naturally related to the final spinor basis.
    using Transformation = typename BasisTransformableTraits<FinalSpinorBasis>::Transformation;

    // The type of Jacobi rotation that is naturally related to the final spinor basis.
    using JacobiRotationType = typename JacobiRotatableTraits<FinalSpinorBasis>::JacobiRotationType;

    // The second-quantized representation of the overlap operator related to the final spinor basis.
    using SQOverlapOperator = typename SpinorBasisTraits<FinalSpinorBasis>::SQOverlapOperator;


protected:
    // The transformation that relates the current set of spinors with the atomic spinors.
    Transformation C;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param C                The transformation that relates the current set of spinors with the atomic spinors.
     */
    SimpleSpinorBasis(const Transformation& C) :
        C {C} {}


    /* 
     *  MARK: Coefficient access
     */

    /**
     *  @return A read-only reference to the transformation that relates the current set of spinors with the atomic spinors.
     */
    const Transformation& expansion() const { return this->C; }


    /*
     *  MARK: General information
     */

    /**
     *  @return The dimension of this simple spinor basis, i.e. the dimension of the underlying coefficient matrix.
     */
    size_t simpleDimension() const { return this->expansion().dimension(); }


    /*
     *  MARK: Orthonormality
     */

    /**
     *  @return The overlap (one-electron) operator of this spinor basis.
     */
    SQOverlapOperator overlap() const { return this->derived().quantize(Operator::Overlap()); }

    /**
     *  Check if this spinor basis is orthonormal within the given precision.
     * 
     *  @param precision                The precision used to test orthonormality.
     * 
     *  @return If this spinor basis is orthonormal.
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {

        const auto S = this->overlap().parameters();

        const auto dim = this->simpleDimension();
        return S.isApprox(SquareMatrix<ExpansionScalar>::Identity(dim), precision);
    }

    /**
     *  @return The transformation to the Löwdin basis: T = S_current^{-1/2}.
     */
    Transformation lowdinOrthonormalization() const {

        // Calculate S^{-1/2}, where S is expressed with respect to the current spinors.
        const auto S = this->overlap().parameters();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver {S};
        return Transformation {eigensolver.operatorInverseSqrt()};
    }

    /**
     *  Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix.
     */
    void lowdinOrthonormalize() { this->C = this->lowdinOrthonormalization(); }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The basis-transformed object.
     */
    FinalSpinorBasis transformed(const Transformation& T) const override {

        auto result = this->derived();
        result.C.transform(T);
        return result;
    }

    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<FinalSpinorBasis>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<FinalSpinorBasis>::rotated;


    /*
     *  MARK: Conforming to `JacobiRotatable`.
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The Jacobi-rotated object.
     */
    FinalSpinorBasis rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto J = Transformation::FromJacobi(jacobi_rotation, this->simpleDimension());
        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<FinalSpinorBasis>::rotate;
};


}  // namespace GQCP
