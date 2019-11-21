// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#pragma once


#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"


namespace GQCP {


/**
 *  A base class that represents a spinor basis with alpha and beta components in separate (possibly different) scalar bases
 * 
 *  @tparam _Shell                  the type of shell the underlying scalar bases contain
 */
template <typename _Shell>
class CompoundSpinorBasis {
public:
    using Shell = _Shell;


protected:
    std::array<ScalarBasis<Shell>, 2> scalar_bases;  // the scalar bases for the alpha and beta components


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param alpha_scalar_basis           the scalar basis in which the alpha components are expanded
     *  @param beta_scalar_basis            the scalar basis in which the beta components are expanded
     */
    CompoundSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        scalar_bases ({alpha_scalar_basis, beta_scalar_basis})
    {}


    /**
     *  Construct a compound spinor basis in which both underlying scalar bases are equal
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    CompoundSpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        CompoundSpinorBasis(scalar_basis, scalar_basis)
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param component        the spin component
     * 
     *  @return the scalar basis in which the requested components are expanded
     */
    const ScalarBasis<Shell>& scalarBasis(const SpinComponent& component) const { return this->scalar_bases[component]; }

    /**
     *  @param component        the spin component
     * 
     *  @return the scalar basis in which the requested components are expanded
     */
    size_t numberOfCoefficients(const SpinComponent& component) const { return this->scalarBasis(component).numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this compound spinor basis
     */
    size_t numberOfSpinors() const { 
        
        const auto K_alpha = this->numberOfCoefficients(SpinComponent::ALPHA);
        const auto K_beta = this->numberOfCoefficients(SpinComponent::BETA);

        return K_alpha + K_beta;
    }
};


}  // namespace GQCP
