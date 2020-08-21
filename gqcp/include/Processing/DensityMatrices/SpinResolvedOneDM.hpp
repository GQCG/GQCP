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
#include "Processing/DensityMatrices/OneDM.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha-alpha and beta-beta spin-resolved density matrices.
 *
 *  @tparam _Scalar             the scalar type of one of the elements
 */
template <typename _Scalar>
class SpinResolvedOneDM {
public:
    using Scalar = _Scalar;


private:
    OneDM<Scalar> D_aa;  // the alpha-alpha 1-DM
    OneDM<Scalar> D_bb;  // the beta-beta 1-DM


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Create a SpinResolvedOneDM from its members.
     *
     *  @param D_aa             the alpha-alpha 1-DM
     *  @param D_bb             the beta-beta 1-DM
     */
    SpinResolvedOneDM(const OneDM<Scalar>& D_aa, const OneDM<Scalar>& D_bb) :
        D_aa {D_aa},
        D_bb {D_bb} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create a spin-resolved 1-DM as half of the total 1-DM.
     * 
     *  @param D            the spin-summed 1-DM
     * 
     *  @return a spin-resolved 1-DM
     */
    static SpinResolvedOneDM<Scalar> FromRestricted(const OneDM<Scalar>& D) {

        const auto D_half = D / 2;
        return SpinResolvedOneDM<Scalar>(D_half, D_half);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the alpha-alpha part of the spin-resolved 1-DM
     */
    const OneDM<Scalar>& alpha() const { return this->D_aa; }

    /**
     *  @return the beta-beta part of the spin-resolved 1-DM
     */
    const OneDM<Scalar>& beta() const { return this->D_bb; }

    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of orbitals that correspond to the given spin
     */
    size_t numberOfOrbitals(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->alpha().dimension();
        }
        case Spin::beta: {
            return this->beta().dimension();
        }
        }
    }


    /**
     *  @return the spin-density matrix, i.e. the difference between the alpha and beta 1-DM
     */
    OneDM<Scalar> spinDensity() const {
        return this->alpha() - this->beta();
    }

    /**
     *  @return the spin-summed density matrix, i.e. the sum of the alpha and beta 1-DM
     */
    OneDM<Scalar> spinSummed() const {
        return this->alpha() + this->beta();
    }
};


}  // namespace GQCP
