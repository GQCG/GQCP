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


#include "DensityMatrix/TwoDM.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  A type that encapsulates the parts of the spin-resolved two-electron density matrix.
 *
 *  @tparam _Scalar         the scalar type of one of the matrix elements
 */
template <typename _Scalar>
struct SpinResolvedTwoDM {
public:
    using Scalar = _Scalar;


private:
    TwoDM<Scalar> d_aaaa;  // the alpha-alpha-alpha-alpha 2-DM
    TwoDM<Scalar> d_aabb;  // the alpha-alpha-beta-beta 2-DM
    TwoDM<Scalar> d_bbaa;  // the beta-beta-alpha-alpha 2-DM
    TwoDM<Scalar> d_bbbb;  // the beta-beta-beta-beta 2-DM


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct a spin-resolved 2-DM from its members.
     *
     *  @param d_aaaa       the alpha-alpha-alpha-alpha 2-DM
     *  @param d_aabb       the alpha-alpha-beta-beta 2-DM
     *  @param d_bbaa       the beta-beta-alpha-alpha 2-DM
     *  @param d_bbbb       the beta-beta-beta-beta 2-DM
     */
    SpinResolvedTwoDM(const TwoDM<Scalar>& d_aaaa, const TwoDM<Scalar>& d_aabb, const TwoDM<Scalar>& d_bbaa, const TwoDM<Scalar>& d_bbbb) :
        d_aaaa {d_aaaa},
        d_aabb {d_aabb},
        d_bbaa {d_bbaa},
        d_bbbb {d_bbbb} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the alpha-alpha-alpha-alpha part of the spin-resolved 2-DM
     */
    const TwoDM<Scalar>& alphaAlpha() const { return this->d_aaaa; }

    /**
     *  @return the alpha-alpha-beta-beta part of the spin-resolved 2-DM
     */
    const TwoDM<Scalar>& alphaBeta() const { return this->d_aabb; }

    /**
     *  @return the beta-beta-alpha-alpha part of the spin-resolved 2-DM
     */
    const TwoDM<Scalar>& betaAlpha() const { return this->d_bbaa; }

    /**
     *  @return the beta-beta-beta-beta part of the spin-resolved 2-DM
     */
    const TwoDM<Scalar>& betaBeta() const { return this->d_bbbb; }

    /**
     *  @param sigma            alpha or beta
     *  @param tau              alpha or beta
     * 
     *  @return the number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the spin-resolved 2-DM
     */
    size_t numberOfOrbitals(const Spin sigma, const Spin tau) const {

        if (sigma == Spin::alpha && tau == Spin::alpha) {
            return this->alphaAlpha().numberOfOrbitals();
        } else if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaBeta().numberOfOrbitals();
        } else if (sigma == Spin::beta && tau == Spin::alpha) {
            return this->betaAlpha().numberOfOrbitals();
        } else {
            return this->betaBeta().numberOfOrbitals();
        }
    }


    /**
     *  @return The orbital (total, spin-summed) two-electron density matrix.
     */
    TwoDM<Scalar> orbitalDensity() const {
        return TwoDM<Scalar>(this->alphaAlpha().Eigen() + this->alphaBeta().Eigen() + this->betaAlpha().Eigen() + this->betaBeta().Eigen());
    }
};


}  // namespace GQCP
