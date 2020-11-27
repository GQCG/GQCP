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


#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Basis/Transformations/UJacobiRotation.hpp"


namespace GQCP {


/**
 *  An interface that implements conformance to `JacobiRotatable` for doubly spin-resolved types.
 * 
 *  @tparam T           The spin-resolved type that should conform to `JacobiRotatable`.
 */
template <typename T>
class DoublySpinResolvedJacobiRotatable:
    public JacobiRotatable<T> {

public:
    // The type of Jacobi rotation for which the basis rotation should be defined.
    using JacobiRotationType = UJacobiRotation;


public:
    /*
     *  MARK: Conforming to `JacobiRotatable`
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     * 
     *  @param jacobi_rotation          The Jacobi rotation.
     * 
     *  @return The jacobi-transformed object.
     */
    T rotated(const JacobiRotationType& jacobi_rotation) const override {

        // Transform the components of 'this' with the components of the transformation matrix.
        const auto transformed_aa = static_cast<const T&>(*this).alphaAlpha().rotated(jacobi_rotation.alpha());

        auto transformed_ab = static_cast<const T&>(*this).alphaBeta().rotated(jacobi_rotation.alpha(), Spin::alpha);
        transformed_ab.rotate(jacobi_rotation.beta(), Spin::beta);

        auto transformed_ba = static_cast<const T&>(*this).betaAlpha().rotated(jacobi_rotation.beta(), Spin::alpha);
        transformed_ba.rotate(jacobi_rotation.alpha(), Spin::beta);

        const auto transformed_bb = static_cast<const T&>(*this).betaBeta().rotated(jacobi_rotation.beta());


        return T {transformed_aa, transformed_ab, transformed_ba, transformed_bb};
    }
};


}  // namespace GQCP
