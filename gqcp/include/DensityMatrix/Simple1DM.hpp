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
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/*
 *  MARK: Simple1DM implementation
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
    public BasisTransformable<Simple1DM<_Scalar, _DerivedDM>> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedDM = _DerivedDM;

    // The type of 'this'.
    using Self = Simple1DM<Scalar, DerivedDM>;

    // The type of transformation that is naturally related to the DerivedDM.
    using Transformation = typename DensityMatrixTraits<DerivedDM>::Transformation;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit base constructors.
    using SquareMatrix<Scalar>::SquareMatrix;


    /*
     *  MARK: General information
     */
    size_t numberOfOrbitals() const { return this->dimension(); }


    /*
     *  MARK: Transformations
     */

    /**
     *  Apply the basis transformation and return the resulting 1-DM.
     * 
     *  @param T        The basis transformation.
     * 
     *  @return The basis-transformed 1-DM.
     */
    Self transformed(const Transformation& T) const override {

        return Self(T.matrix().inverse().conjugate() * (*this) * T.matrix().inverse().transpose());  // Note that this basis transformation formula is different from the one-electron operator one. See `SimpleSQOneElectronOperator`.
    }
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename Scalar, typename DerivedDM>
struct BasisTransformableTraits<Simple1DM<Scalar, DerivedDM>> {

    // The type of the transformation for which the basis transformation should be defined.
    using Transformation = typename DensityMatrixTraits<DerivedDM>::Transformation;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename Scalar, typename DerivedDM>
struct JacobiRotatableTraits<Simple1DM<Scalar, DerivedDM>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
