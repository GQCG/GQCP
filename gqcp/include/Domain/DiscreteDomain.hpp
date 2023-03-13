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
#include "Domain/SimpleDomain.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/CRTP.hpp"
#include "Utilities/type_traits.hpp"

#include <boost/dynamic_bitset.hpp>

#include <algorithm>


namespace GQCP {


/**
 *
 */
class DiscreteDomain:
    public SimpleDomain<DiscreteDomain> {
private:
    // The representation of the discrete domain as a bitstring.
    boost::dynamic_bitset<> domain;
    // A set of indices that correspond to the elements included in the discrete domain.
    std::vector<size_t> domain_indices;

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a discrete domain from a set of domain indices.
     *
     *  @param domain_indices             The indices that the domain contains.
     *  @param M                            The dimension of the discrete domain, i.e. the maximum number of elements that the discrete domain can contain.
     *
     *  @return A spin-unresolved ONV from a set of occupied indices.
     */
    DiscreteDomain(const std::vector<size_t>& domain_indices, const size_t M);

    /**
     *  Create a discrete domain from an unsigned representation.
     *
     *  @param unsigned_representation          The representation of this discrete domain as an unsigned integer.
     *  @param M                                The dimension of the discrete domain, i.e. the maximum number of elements that the discrete domain can contain.
     */
    DiscreteDomain(const size_t unsigned_representation, const size_t M);

    /**
     *  Add an element to the domain at position i.
     *
     *  @param i        The index in the domain bitstring.
     */
    void addElement(const size_t i);

    /**
     *  @return The domain representation as a bitstring.
     */
    const boost::dynamic_bitset<>& asBitstring() const { return this->domain; }

    /**
     *  @return The domain string as a bitstring.
     */
    std::string asString() const;

    /**
     *  @return     The dimension of this discrete domain, i.e. the maximum number of elements that the discrete domain can contain.
     */
    size_t dimension() const { return this->domain.size(); }


    const std::vector<size_t>& domainIndices() const { return this->domain_indices; }

    /**
     *  @param i        The index in the domain bitstring.
     *
     *  @return whether the Domain contains an element at index i.
     */
    bool inDomain(const size_t& i) const { return this->domain[i]; }

    /**
     *  @return     The number of elements the domain contains.
     */
    size_t numberOfElements() const { return this->domain_indices.size(); }

    /**
     *  @param rhs        The other CartesianGTO.
     *
     *  @return     If this discrete domain is equal to the other discrete domain.
     */
    bool operator==(const DiscreteDomain& rhs) const;

    /**
     *  @param i            The domain index.
     *
     *  @return     The i-th element with values 0 or 1 depending on whether the i-th element is in the domain.
     */
    bool operator()(const size_t i) const { return this->domain[i]; }

    /**
     *  @param rhs            The discrete domain that will be assigned to `this'.
     *
     *  @return     A discrete domain with the class variables provided by `rhs'.
     */
    DiscreteDomain& operator=(const DiscreteDomain& rhs);


    /**
     * Calculate the overlap between two discrete domains, i.e. the number of matching domain elements.
     *
     *  @param other            The other discrete domain.
     *
     *  @return     The number of domain elements that are equal between `this' and the other discrete domain.
     */
    size_t overlapWith(const DiscreteDomain& other) const;

    /**
     *  Remove an element from the domain at position i.
     *
     *  @param i        The index in the domain bitstring.
     */
    void removeElement(const size_t i);

    /**
     *  @return The unsigned representation of this domain.
     */
    size_t unsignedRepresentation() const { return this->domain.to_ulong(); }
};


}  // namespace GQCP
