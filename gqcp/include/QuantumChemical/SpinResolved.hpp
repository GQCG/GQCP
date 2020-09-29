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


#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  A utility type encapsulating an alpha- and beta-type.
 * 
 *  @param _Of          The type of the alpha- and beta-objects. The name 'Of' is chosen for a natural reading `SpinResolved<_Of, _Derived>`.
 *  @param _Derived     The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
 */
template <typename _Of, typename _Derived>
class SpinResolved {
public:
    // The type of the alpha- and beta-objects. The name 'Of' is chosen for a natural reading `SpinResolved<_Of, _Derived>`.
    using Of = _Of;

    // The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
    using Derived = _Derived;

    // The type of 'this'.
    using Self = SpinResolved<Of, Derived>;


private:
    // The alpha-object.
    Of m_alpha;

    // The beta-object.
    Of m_beta;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a spin-resolved type from its alpha and beta objects.
     * 
     *  @param alpha        The alpha-object.
     *  @param beta         The beta-object.
     */
    SpinResolved(const Of& alpha, const Of& beta) :
        m_alpha {alpha},
        m_beta {beta} {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a derived (spin-resolved) type, from equal alpha and beta representations.
     * 
     *  @param equal        The equal representation for both the alpha and beta objects.
     * 
     *  @return The dervived (spin-resolved) type.
     */
    static Derived FromEqual(const Of& equal) {
        return Derived {equal, equal};
    }


    /*
     *  MARK: Accessing spin components
     */

    /**
     *  @return A read-only reference to the alpha object.
     */
    const Of& alpha() const { return this->m_alpha; }

    /**
     *  @return A writable reference to the alpha object.
     */
    Of& alpha() { return this->m_alpha; }

    /**
     *  @return A read-only reference to the beta object.
     */
    const Of& beta() const { return this->m_beta; }

    /**
     *  @return A writable reference to the beta object.
     */
    Of& beta() { return this->m_beta; }

    /**
     *  Access the alpha or beta component of this spin-resolved type.
     * 
     *  @param sigma            Alpha or beta.
     * 
     *  @return A read-only reference to the alpha or beta object.
     */
    const Of& component(const Spin sigma) const {

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

    /**
     *  Access the alpha or beta component of this spin-resolved type.
     * 
     *  @param sigma            Alpha or beta.
     * 
     *  @return A writable reference to the alpha or beta object.
     */
    Of& component(const Spin sigma) { return const_cast<Of&>(const_cast<const Self*>(this)->component(sigma)); }
};


}  // namespace GQCP
