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
#include "Basis/Transformations/JacobiRotation.hpp"
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/*
 *  MARK: `Simple1DM` implementation
 */

/**
 *  A one-electron density matrix that is described by a single matrix.
 * 
 *  This class is used as a base class for `Orbital1DM` and `G1DM`, since they are both expressed using a single matrix, as opposed to `SpinResolved1DM`, which uses separate alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                 The scalar type used for a density matrix element: real or complex.
 *  @tparam _DerivedDM              The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedDM>
class Simple1DM:
    public SquareMatrix<_Scalar>,
    public BasisTransformable<_DerivedDM>,
    public JacobiRotatable<_DerivedDM> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedDM = _DerivedDM;

    // The type of 'this'.
    using Self = Simple1DM<Scalar, DerivedDM>;

    // The type of transformation that is naturally related to the `DerivedDM`.
    using Transformation = typename DensityMatrixTraits<DerivedDM>::Transformation;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SquareMatrix`' constructors.
    using SquareMatrix<Scalar>::SquareMatrix;


    /*
     *  MARK: General information
     */
    size_t numberOfOrbitals() const { return this->dimension(); }


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the resulting 1-DM.
     * 
     *  @param T        The basis transformation.
     * 
     *  @return The basis-transformed 1-DM.
     */
    DerivedDM transformed(const Transformation& T) const override {

        // The transformation formulas for one-electron operators and 1-DMs are similar, but not quite the same. Instead of using T, the transformation formula for the 1-DM uses T_inverse_transpose. See also (https://gqcg-res.github.io/knowdes/spinor-transformations.html).
        const auto T_related = T.matrix().transpose().inverse();
        return DerivedDM {T_related.adjoint() * (*this) * T_related};
    }

    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedDM>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedDM>::rotated;


    /*
     *  MARK: Conforming to `JacobiRotatable`
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The Jacobi-transformed object.
     */
    DerivedDM rotated(const JacobiRotation& jacobi_rotation) const override {

        // We implement this rotation by constructing a Jacobi rotation matrix and then simply doing a rotation with it.
        const auto J = Transformation::FromJacobi(jacobi_rotation, this->numberOfOrbitals());
        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedDM>::rotate;
};


}  // namespace GQCP
