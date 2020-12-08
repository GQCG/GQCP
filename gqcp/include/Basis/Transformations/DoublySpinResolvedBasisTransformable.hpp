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
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  An interface that implements conformance to `BasisTransformable` for doubly spin-resolved types.
 * 
 *  @tparam Type            The spin-resolved type that should conform to `BasisTransformable`.
 */
template <typename Type>
class DoublySpinResolvedBasisTransformable:
    public BasisTransformable<Type> {

public:
    // The type of the transformation for which the basis transformation should be defined.
    using Transformation = typename BasisTransformableTraits<Type>::Transformation;


public:
    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     * 
     *  @param T            The basis transformation.
     * 
     *  @return The basis-transformed spin-resolved object.
     */
    Type transformed(const Transformation& T) const override {

        // Transform the components of 'this' with the components of the transformation matrix.
        const auto transformed_aa = static_cast<const Type&>(*this).alphaAlpha().transformed(T.alpha());

        auto transformed_ab = static_cast<const Type&>(*this).alphaBeta().transformed(T.alpha(), Spin::alpha);
        transformed_ab.transform(T.beta(), Spin::beta);

        auto transformed_ba = static_cast<const Type&>(*this).betaAlpha().transformed(T.beta(), Spin::alpha);
        transformed_ba.transform(T.alpha(), Spin::beta);

        const auto transformed_bb = static_cast<const Type&>(*this).betaBeta().transformed(T.beta());


        return Type {transformed_aa, transformed_ab, transformed_ba, transformed_bb};
    }
};


}  // namespace GQCP
