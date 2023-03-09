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


#include "Domain/DiscreteDomain.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *
 *
 *
 */
DiscreteDomain::DiscreteDomain(const VectorX<size_t>& indices, const size_t M) :
    m_indices {indices} {
    // Generate the corresponding unsigned representation associated with the indices of the domain.
    size_t unsigned_representation = 0;
    for (size_t i = 0; ++i; i < m_indices.size()) {
        unsigned_representation += std::pow(2, m_indices[i]);
    }
    this->domain = boost::dynamic_bitset<> {M, unsigned_representation};
}


/**
 *  Add an element to the domain at position i.
 *
 *  @param i        The index in the domain bitstring.
 */
void DiscreteDomain::addElement(const size_t i) {
    if (i < 0) {
        throw std::invalid_argument("DiscreteDomain::addElement(const size_t): Domain-index is smaller than zero. Please provide a valid index.");
    }
    if (i >= this->dimension()) {
        throw std::invalid_argument("DiscreteDomain::addElement(const size_t): Domain-index is larger than the domain size. Provide a valid index.");
    }
    if (this->domain[i] == 1) {  // Check whether the element at index i is already present in the domain.
        throw std::invalid_argument("DiscreteDomain::addElement(const size_t): You cannot add an element at index i to the discrete domain if the element is already present.");
    }

    this->domain.set(i, true);
}


/**
 *  @return The domain string as a bitstring.
 */
std::string DiscreteDomain::asString() const {
    const auto K = this->dimension();

    std::string s;
    for (size_t i = 0; i < K; ++i) {
        s += std::to_string(this->operator[](i));
    }
    return s;
}


/**
 *  @param other        The other CartesianGTO.
 *
 *  @return if this discrete domain is equal to the other discrete domain.
 */
bool DiscreteDomain::operator==(const DiscreteDomain& other) const {
    return boost::operator==(this->domain, other.asBitstring());
}


/**
 *  Remove an element from the domain at position i.
 *
 *  @param i        The index in the domain bitstring.
 */
void DiscreteDomain::removeElement(const size_t i) {
    if (i < 0) {
        throw std::invalid_argument("DiscreteDomain::removeElement(const size_t): Domain-index is smaller than zero. Please provide a valid index.");
    }
    if (i >= this->dimension()) {
        throw std::invalid_argument("DiscreteDomain::removeElement(const size_t): Domain-index is larger than the domain size. Provide a valid index.");
    }
    if (this->domain[i] == 0) {  // Check whether the element at index i is actually present in the domain.
        throw std::invalid_argument("DiscreteDomain::removeElement(const size_t): You cannot remove an element at index i from the discrete domain if there is no element at that index.");
    }

    this->domain.set(i, false);
}


/*
 *  MARK: DomainTraits
 */


}  // namespace GQCP
