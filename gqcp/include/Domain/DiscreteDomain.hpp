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


#include "Domain/SimpleDomain.hpp"
#include "Utilities/CRTP.hpp"

#include <boost/dynamic_bitset.hpp>


namespace GQCP {


/**
 *  A type specifically designed to act as a parent class for `DiscreteDomain` and `ContinuousDomain` in order to share common functionality.
 */
class DiscreteDomain:
    public SimpleDomain<DiscreteDomain> {
private:
    // The representation of the discrete domain as a bitstring.
    boost::dynamic_bitset<> domain;
    // A set of indices that correspond to the elements included in the discrete domain.
    VectorX<size_t> m_indices;

public:
    DiscreteDomain(const std::vector<size_t>& indices, const size_t M) :
        m_indices {indices} {
        // Generate the corresponding unsigned representation associated with the indices of the domain.
        size_t unsigned_representation = 0;
        for (const auto& index : indices) {
            unsigned_representation += std::pow(2, index);
        }
        this->domain = boost::dynamic_bitset<> { M, unsigned_representation }
    }

    /**
     *  @return the domain representation as a bitstring.
     */
    const boost::dynamic_bitset<>& asBitstring() const { return this->domain; }

    /**
     *  @param element        the element in the domain
     *
     *  @return whether the Domain contains `element'.
     */
    bool inDomain(const size_t element) const { return this->domain[element]; }

    /**
     *  @param other        the other CartesianGTO
     *
     *  @return if this DerivedDomain is equal to the other DerivedDomain.
     */
    virtual bool operator==(const DiscreteDomain& other) const { return boost::operator==(this->domain, other.asBitstring()); }
};


/*
 *  MARK: DomainTraits
 */

/**
 *  A type that provides compile-time information on `DiscreteDomain` that is otherwise not accessible through a public class alias.
 */
struct DomainTraits<DiscreteDomain> {

    // The type of elements that are present in the domain.
    using ElementType = size_t;
};


}  // namespace GQCP
