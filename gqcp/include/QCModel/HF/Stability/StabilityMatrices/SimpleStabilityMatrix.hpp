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


#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/*
 *  MARK: SimpleStabilityMatrix implementation
 */

/**
 *  A stability matrix of a Hartree-Fock QC model
 * 
 *  This class is used as a base class for the different stability matrices of the Hartree-Fock methods. The word 'simple' is used here as an antonym for 'compound'.
 * 
 *  @tparam _Scalar                              The scalar type used for a stability matrix element: real or complex.
 *  @tparam _DerivedStabilityMatrix              The type of the stability matrix that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _Scalar, typename _DerivedStabilityMatrix>
class SimpleStabilityMatrix:
    public Matrix<_Scalar> {
public:
    // The scalar type used for a density matrix element: real or complex.
    using Scalar = _Scalar;

    // The type of the density matrix that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedStabilityMatrix = _DerivedStabilityMatrix;

    // The type of 'this'.
    using Self = SimpleStabilityMatrix<Scalar, DerivedStabilityMatrix>;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit base constructors.
    using Matrix<Scalar>::Matrix;
};

}  // namespace GQCP
