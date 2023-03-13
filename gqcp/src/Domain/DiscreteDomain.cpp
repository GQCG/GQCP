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
 *  Create a discrete domain from a set of domain indices.
 *
 *  @param indices                      The indices that the domain contains.
 *  @param M                            The dimension of the discrete domain, i.e. the maximum number of elements that the discrete domain can contain.
 *
 *  @return A spin-unresolved ONV from a set of occupied indices.
 */
DiscreteDomain::DiscreteDomain(const std::vector<size_t>& indices, const size_t M) :
    domain_indices {indices} {

    // Generate the corresponding unsigned representation associated with the indices of the domain.
    size_t unsigned_representation = 0;
    for (size_t i = 0; i < domain_indices.size(); i++) {
        unsigned_representation += std::pow(2, domain_indices[i]);
    }
    this->domain = boost::dynamic_bitset<> {M, unsigned_representation};
}


/**
 *  Create a discrete domain from an unsigned representation.
 *
 *  @param unsigned_representation          The representation of this discrete domain as an unsigned integer.
 *  @param M                                The dimension of the discrete domain, i.e. the maximum number of elements that the discrete domain can contain.
 */
DiscreteDomain::DiscreteDomain(size_t unsigned_representation, const size_t M) :
    domain {boost::dynamic_bitset<> {M, unsigned_representation}} {

    std::vector<size_t> domain_indices {};
    while (unsigned_representation != 0) {
        domain_indices.push_back(__builtin_ctzl(unsigned_representation));                // Retrieves the domain index. Note: bitstring is read from right to left, domain index from left to right.
        unsigned_representation ^= (unsigned_representation & -unsigned_representation);  // Negative of an unsigned type: http://www.cplusplus.com/forum/beginner/182082/.
    }

    this->domain_indices = domain_indices;
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
    this->domain_indices.push_back(i);
    sort(this->domain_indices.begin(), this->domain_indices.end());
}


/**
 *  @return The domain string as a bitstring.
 */
std::string DiscreteDomain::asString() const {
    std::string text;  // The string that will contain the textual representation of this discrete domain.

    boost::to_string(this->domain, text);
    std::reverse(text.begin(), text.end());  // The bitstring is read from right to left, but we return it from left to right for easier reading.

    return text;
}


/**
 *  @param other        The other CartesianGTO.
 *
 *  @return if this discrete domain is equal to the other discrete domain.
 */
bool DiscreteDomain::operator==(const DiscreteDomain& rhs) const {
    return boost::operator==(this->domain, rhs.asBitstring());
}


/**
 *  @param rhs            The discrete domain that will be assigned to `this`.
 *
 *  @return     A discrete domain with the class variables provided by `rhs'.
 */
DiscreteDomain& DiscreteDomain::operator=(const DiscreteDomain& rhs) {
    this->domain = rhs.asBitstring();
    this->domain_indices = rhs.domainIndices();

    return (*this);
}


/**
 *  @param os       The output stream which the discrete domain should be concatenated to.
 *  @param domain      The discrete domain that should be concatenated to the output stream.
 *
 *  @return The updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const DiscreteDomain& domain) {
    return os << domain.asString();
}


/**
 * Calculate the overlap between two discrete domains, i.e. the number of matching domain elements.
 *
 *  @param other            The other discrete domain.
 *
 *  @return     The number of domain elements that are equal between `this` and the other discrete domain.
 */
size_t DiscreteDomain::overlapWith(const DiscreteDomain& other) const {

    boost::dynamic_bitset<> overlap;
    overlap = this->domain & other.asBitstring();

    return overlap.count();
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
    this->domain_indices.erase(std::remove(this->domain_indices.begin(), this->domain_indices.end(), i), this->domain_indices.end());
}


}  // namespace GQCP
