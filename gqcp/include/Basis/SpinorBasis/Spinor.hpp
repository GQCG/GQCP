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
#include "Mathematical/AbstractFunction/LinearCombination.hpp"


namespace GQCP {


/**
 *  @tparam _Scalar                 the scalar representation of the expansion coefficients
 *  @tparam _BasisFunction          the type of the underlying basis functions
 */
template <typename _Scalar, typename _BasisFunction>
class Spinor {
public:
    using Scalar = _Scalar;
    using BasisFunction = _BasisFunction;


private:
    LinearCombination<Scalar, BasisFunction> alpha_component;  // the alpha-component of the spinor, as an expansion of underlying scalar functions
    LinearCombination<Scalar, BasisFunction> beta_component;   // the beta-component of the spinor, as an expansion of underlying scalar functions


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A memberwise constructor.
     * 
     *  @param alpha_component              the alpha-component of the spinor, as an expansion of underlying scalar functions
     *  @param beta_component               the beta-component of the spinor, as an expansion of underlying scalar functions
     */
    Spinor(const LinearCombination<Scalar, BasisFunction>& alpha_component, const LinearCombination<Scalar, BasisFunction>& beta_component) :
        alpha_component {alpha_component},
        beta_component {beta_component} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the sigma-component of this spinor
     */
    const LinearCombination<Scalar, BasisFunction>& component(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha:
            return this->alpha_component;
            break;

        case Spin::beta:
            return this->beta_component;
            break;
        }
    }
};


}  // namespace GQCP
