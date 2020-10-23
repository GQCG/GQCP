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

#include <array>
#include <vector>


namespace GQCP {


/**
 *  A utility type encapsulating four objects, each for every combination of alpha and beta. In contrast to `DoubleSpinResolved`, this type is to be used as a base class/interface for types that are doubly spin-resolved.
 * 
 *  @param _Of          The type of the doubly spin-resolved objects. The name 'Of' is chosen for a natural reading `DoublySpinResolvedBase<_Of, _Derived>`.
 *  @param _Derived     The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
 */
template <typename _Of, typename _Derived>
class DoublySpinResolvedBase {
public:
    // The type of the doubly spin-resolved objects. The name 'Of' is chosen for a natural reading `DoublySpinResolvedBase<_Of, _Derived>`.
    using Of = _Of;

    // The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
    using Derived = _Derived;

    // The type of 'this'.
    using Self = DoublySpinResolvedBase<Of, Derived>;


private:
    // The alpha-alpha-object.
    Of aa;

    // The alpha-beta object.
    Of ab;

    // The beta-alpha object.
    Of ba;

    // The beta-beta object.
    Of bb;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a doubly spin-resolved instance from its constituent objects.
     * 
     *  @param aa       The alpha-alpha-object.
     *  @param ab       The alpha-beta object.
     *  @param ba       The beta-alpha object.
     *  @param bb       The beta-beta object.
     */
    DoublySpinResolvedBase(const Of& aa, const Of& ab, const Of& ba, const Of& bb) :
        aa {aa},
        ab {ab},
        ba {ba},
        bb {bb} {}


    /**
     *  Construct a doubly spin-resolved instance from a vector containing its alpha and beta objects.
     * 
     *  @param all          A vector containing all alpha and beta objects, in the canonical order.
     */
    DoublySpinResolvedBase(const std::vector<Of>& all) :
        DoublySpinResolvedBase(all[0], all[1], all[2], all[3]) {

        if (all.size() != 4) {
            throw std::invalid_argument("DoublySpinResolvedBase(const std::vector<Of>&): The given vector does not have exactly four elements.");
        }
    }


    /**
     *   Construct a doubly spin-resolved instance from an array containing its alpha and beta objects.
     * 
     *  @param all          An array containing all alpha and beta objects, in the canonical order.
     */
    DoublySpinResolvedBase(const std::array<Of, 4>& all) :
        DoublySpinResolvedBase(all[0], all[1], all[2], all[3]) {}


    /**
     *  Construct a doubly spin-resolved type from an initializer list containing its alpha and beta objects.
     * 
     *  @param all          An initializer list containing both all alpha and beta objects (in that order).
     */
    DoublySpinResolvedBase(const std::initializer_list<Of>& all) :
        DoublySpinResolvedBase(std::vector<Of>(all)) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a derived (doubly spin-resolved) type, from equal alpha and beta representations.
     * 
     *  @param equal        The equal representation for all the alpha and beta objects.
     * 
     *  @return The derived (spin-resolved) type.
     */
    static Derived FromEqual(const Of& equal) {
        return Derived {equal, equal, equal, equal};
    }


    /*
     *  MARK: Accessing spin components
     */

    /**
     *  @return A read-only reference to the alpha-alpha object.
     */
    const Of& alphaAlpha() const { return this->aa; }

    /**
     *  @return A writable reference to the alpha-alpha object.
     */
    Of& alphaAlpha() { return this->aa; }

    /**
     *  @return A read-only reference to the alpha-beta object.
     */
    const Of& alphaBeta() const { return this->ab; }

    /**
     *  @return A writable reference to the alpha-beta object.
     */
    Of& alphaBeta() { return this->ab; }

    /**
     *  @return A read-only reference to the beta-alpha object.
     */
    const Of& betaAlpha() const { return this->ba; }

    /**
     *  @return A writable reference to the beta-alpha object.
     */
    Of& betaAlpha() { return this->ba; }

    /**
     *  @return A read-only reference to the beta-beta object.
     */
    const Of& betaBeta() const { return this->bb; }

    /**
     *  @return A writable reference to the beta-beta object.
     */
    Of& betaBeta() { return this->bb; }

    /**
     *  Access one of the components of this doubly spin-resolved instance.
     * 
     *  @param sigma            Alpha or beta.
     *  @param tau              Alpha or beta.
     * 
     *  @return A read-only reference to one of the components.
     */
    const Of& component(const Spin sigma, const Spin tau) const {

        if (sigma == Spin::alpha && tau == Spin::alpha) {
            return this->alphaAlpha();
        } else if (sigma == Spin::alpha && tau == Spin::beta) {
            return this->alphaBeta();
        } else if (sigma == Spin::beta && tau == Spin::alpha) {
            return this->betaAlpha();
        } else {
            return this->betaBeta();
        }
    }

    /**
     *  Access one of the components of this doubly spin-resolved instance.
     * 
     *  @param sigma            Alpha or beta.
     *  @param tau              Alpha or beta.
     * 
     *  @return A writable reference to one of the components.
     */
    Of& component(const Spin sigma, const Spin tau) { return const_cast<Of&>(const_cast<const Self*>(this)->component(sigma, tau)); }
};


}  // namespace GQCP
