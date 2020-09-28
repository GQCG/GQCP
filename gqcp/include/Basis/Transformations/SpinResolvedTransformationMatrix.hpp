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


#include "Basis/SpinorBasis/Spin.hpp"
#include "Basis/Transformations/TransformationMatrix.hpp"


namespace GQCP {


/**
 *  A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases.
 */
template <typename _Scalar>
class SpinResolvedTransformationMatrix {
public:
    using Scalar = _Scalar;


private:
    TransformationMatrix<Scalar> T_a;  // the transformation matrix for the alpha spin-orbitals
    TransformationMatrix<Scalar> T_b;  // the transformation matrix for the beta spin-orbitals

public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Create a spin-resolved transformation matrix from its members.
     * 
     *  @param T_a          the transformation matrix for the alpha spin-orbitals
     *  @param T_b          the transformation matrix for the beta spin-orbitals
     */
    SpinResolvedTransformationMatrix(const TransformationMatrix<Scalar>& T_a, const TransformationMatrix<Scalar>& T_b) :
        T_a {T_a},
        T_b {T_b} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param T                the equal transformation matrix for both the alpha and the beta spin-orbitals
     * 
     *  @return a spin-resolved transformation matrix where the transformation matrices for the alpha and beta spin-orbitals are equal
     */
    static SpinResolvedTransformationMatrix<Scalar> FromRestricted(const TransformationMatrix<Scalar>& T) {

        return SpinResolvedTransformationMatrix<Scalar>(T, T);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the transformation matrix for the alpha spin-orbitals
     */
    const TransformationMatrix<Scalar> alpha() const { return this->T_a; }

    /**
     *  @return the transformation matrix for the beta spin-orbitals
     */
    const TransformationMatrix<Scalar> beta() const { return this->T_b; }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of orbitals (spin-orbitals) that the transformation matrix for the sigma spin-orbitals is related to
     */
    size_t numberOfOrbitals(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->alpha().numberOfOrbitals();
            break;
        }

        case Spin::beta: {
            return this->beta().numberOfOrbitals();
            break;
        }
        }
    }
};


}  // namespace GQCP
