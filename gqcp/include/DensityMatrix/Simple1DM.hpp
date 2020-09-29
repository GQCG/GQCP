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
#include "Mathematical/Representation/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A one-electron density matrix that is described by a single matrix.
 * 
 *  This class is used as a base class for `Orbital1DM` and `G1DM`, since they are both expressed using a single matrix, as opposed to `U1DM`, which uses separate alpha- and beta- matrices. The word 'simple' is used here as an antonym for 'compound'.
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
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed 1-DM.
     */
    Self transformed(const TM& transformation_matrix) const override {

        return Self(transformation_matrix.conjugate() * (*this) * transformation_matrix.transpose());  // Note that this basis transformation formula is different from the one-electron operator one. See SimpleSQOneElectronOperator.
    }
};


}  // namespace GQCP
