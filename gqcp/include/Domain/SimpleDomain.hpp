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


#include "Domain/DomainTraits.hpp"
#include "Utilities/CRTP.hpp"


namespace GQCP {


/**
 *  A type specifically designed to act as a parent class for `DiscreteDomain` and `ContinuousDomain` in order to share common functionality.
 *
 *  @tparam _DerivedDomain              The type of the domain that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _DerivedDomain>
class SimpleDomain:
    public CRTP<_DerivedDomain> {
public:
    // The type of domain that derives from this class, enabling CRTP and compile-time polymorphism..
    using DerivedDomain = _DerivedDomain;
    // The type of elements that are present in the domain: discrete or continuous.
    using ElementType = typename DomainTraits<DerivedDomain>::ElementType;


public:
    /**
     *  @param element        the element in the domain
     *
     *  @return whether the Domain contains `element'.
     */
    virtual bool inDomain(const ElementType& element) const = 0;

    /**
     *  @param other        the other domain
     *
     *  @return if this DerivedDomain is equal to the other DerivedDomain.
     */
    virtual bool operator==(const DerivedDomain& other) const = 0;

    /**
     *  @param other        the other domain
     *
     *  @return if this DerivedDomain is not equal to the other DerivedDomain.
     */
    bool operator!=(DerivedDomain& other) const {
        return !this->operator==(other);
    }
};


}  // namespace GQCP
