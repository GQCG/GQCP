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


namespace GQCP {


/**
 *  An interface that implements conformance to `BasisTransformable` for spin-resolved types.
 * 
 *  @tparam T           The spin-resolved type that should conform to `BasisTransformable`.
 */
template <typename T>
class SpinResolvedBasisTransformable:
    public BasisTransformable<T> {

public:
    // The type of the transformation matrix for which the basis transformation should be defined.
    using TM = typename BasisTransformableTraits<T>::TM;


public:
    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param transformation_matrix        The type that encapsulates the basis transformation coefficients.
     * 
     *  @return The basis-transformed spin-resolved object.
     */
    T transformed(const TM& transformation_matrix) const override {

        // Transform the components of 'this' with the components of the transformation matrix.
        const auto alpha_transformed = static_cast<const T&>(*this).alpha().transformed(transformation_matrix.alpha());
        const auto beta_transformed = static_cast<const T&>(*this).beta().transformed(transformation_matrix.beta());

        return T {alpha_transformed, beta_transformed};
    }
};


}  // namespace GQCP
