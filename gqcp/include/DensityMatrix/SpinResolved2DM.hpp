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


#include "DensityMatrix/MixedSpinResolved2DMComponent.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "DensityMatrix/PureSpinResolved2DMComponent.hpp"
#include "QuantumChemical/DoublySpinResolvedBase.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  A type that encapsulates the spin parts of the spin-resolved two-electron density matrix.
 *
 *  @tparam _Scalar         The scalar type of one of the density matrix elements: real or complex.
 */
template <typename _Scalar>
class SpinResolved2DM:
    public DoublySpinResolvedBase<PureSpinResolved2DMComponent<_Scalar>, MixedSpinResolved2DMComponent<_Scalar>, SpinResolved2DM<_Scalar>> {
public:
    // The scalar type of one of the density matrix elements: real or complex.
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `DoublySpinResolvedBase`'s constructors.
    using DoublySpinResolvedBase<PureSpinResolved2DMComponent<_Scalar>, MixedSpinResolved2DMComponent<_Scalar>, SpinResolved2DM<_Scalar>>::DoublySpinResolvedBase;


    /*
     *  MARK: General information
     */

    /**
     *  @param sigma            Alpha or beta.
     *  @param tau              Alpha or beta.
     * 
     *  @return The number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the spin-resolved 2-DM.
     */
    size_t numberOfOrbitals(const Spin sigma, const Spin tau) const {

        if (sigma == Spin::alpha && tau == Spin::beta) {
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
    Orbital2DM<Scalar> orbitalDensity() const {
        return Orbital2DM<Scalar> {this->alphaAlpha().Eigen() + this->alphaBeta().Eigen() + this->betaAlpha().Eigen() + this->betaBeta().Eigen()};
    }
};


}  // namespace GQCP
