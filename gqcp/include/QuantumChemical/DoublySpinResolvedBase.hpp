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


#include <array>
#include <vector>


namespace GQCP {


/**
 *  A utility type encapsulating four objects, each for every combination of alpha and beta.
 *
 *  @param _Pure        The type that represents a 'pure' combination of spin components, i.e. the alpha-alpha or beta-beta type.
 *  @param _Mixed       The type that represents a 'mixed' combination of spin components, i.e. the alpha-beta or beta-alpha type.
 *  @param _Derived     The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
 */
template <typename _Pure, typename _Mixed, typename _Derived>
class DoublySpinResolvedBase {
public:
    // The type that represents a 'pure' combination of spin components, i.e. the alpha-alpha or beta-beta type.
    using Pure = _Pure;

    // The type that represents a 'mixed' combination of spin components, i.e. the alpha-beta or beta-alpha type.
    using Mixed = _Mixed;

    // The type that derives from this type, given as a template argument, enabling CRTP and compile-time polymorphism.
    using Derived = _Derived;

    // The type of 'this'.
    using Self = DoublySpinResolvedBase<Pure, Mixed, Derived>;


private:
    // The alpha-alpha-object.
    Pure aa;

    // The alpha-beta object.
    Mixed ab;

    // The beta-alpha object.
    Mixed ba;

    // The beta-beta object.
    Pure bb;


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
    DoublySpinResolvedBase(const Pure& aa, const Mixed& ab, const Mixed& ba, const Pure& bb) :
        aa {aa},
        ab {ab},
        ba {ba},
        bb {bb} {}


    /*
     *  MARK: Accessing spin components
     */

    /**
     *  @return A read-only reference to the alpha-alpha object.
     */
    const Pure& alphaAlpha() const { return this->aa; }

    /**
     *  @return A writable reference to the alpha-alpha object.
     */
    Pure& alphaAlpha() { return this->aa; }

    /**
     *  @return A read-only reference to the alpha-beta object.
     */
    const Mixed& alphaBeta() const { return this->ab; }

    /**
     *  @return A writable reference to the alpha-beta object.
     */
    Mixed& alphaBeta() { return this->ab; }

    /**
     *  @return A read-only reference to the beta-alpha object.
     */
    const Mixed& betaAlpha() const { return this->ba; }

    /**
     *  @return A writable reference to the beta-alpha object.
     */
    Mixed& betaAlpha() { return this->ba; }

    /**
     *  @return A read-only reference to the beta-beta object.
     */
    const Pure& betaBeta() const { return this->bb; }

    /**
     *  @return A writable reference to the beta-beta object.
     */
    Pure& betaBeta() { return this->bb; }

    /**
     *  @param sigma     The spin sigma for which the pure component is asked.
     *
     *  @return A read-only reference to the pure alpha-alpha or beta-beta object.
     */
    const Pure& pureComponent(const Spin& sigma) const {

        if (sigma == Spin::alpha) {
            return this->aa;
        } else {
            return this->bb;
        }
    }

    /**
     *  @param sigma     The spin sigma for which the pure component is asked.
     *
     *  @return A writable reference to the pure alpha-alpha or beta-beta object.
     */
    Pure& pureComponent(const Spin& sigma) {

        if (sigma == Spin::alpha) {
            return this->aa;
        } else {
            return this->bb;
        }
    }
};


}  // namespace GQCP
