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


#include "Partition/DiscreteDomainPartition.hpp"


namespace GQCP {


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


#include "Domain/DiscreteDomain.hpp"
#include "Partition/DomainPartition.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Create a discrete domain partition from a vector of discrete domains.
 *
 *  @param domains          The vector of discrete domains.
 */
DiscreteDomainPartition::DiscreteDomainPartition(const std::vector<DiscreteDomain>& domains) :
    partition {domains} {
    size_t domain_occupations = 0;

    for (size_t i = 0; i < domains.size() - 1; i++) {
        for (size_t j = i + 1; j < domains.size(); j++) {
            if (domains[i].dimension() != domains[j].dimension()) {
                throw std::invalid_argument("DiscreteDomainPartition::DiscreteDomainPartition(const std::vector<DiscreteDomain> domains): Found domains where the domain dimensions are not of equal size.");
            }
            if (domains[i].overlapWith(domains[j])) {
                throw std::invalid_argument("DiscreteDomainPartition::DiscreteDomainPartition(const std::vector<DiscreteDomain> domains): Found fuzzy domains, please provide domains that do not overlap with eachother.");
            }
        }
        domain_occupations += domains[i].numberOfElements();
    }
    domain_occupations += domains.back().numberOfElements();

    if (domain_occupations != domains[0].dimension()) {
        throw std::invalid_argument("DiscreteDomainPartition::DiscreteDomainPartition(const std::vector<DiscreteDomain> domains): Found non-complete domains, please provide a collection of domains where each index belongs to a domain once.");
    }
}

/*
 *  MARK: General info
 */

/**
 *  @return The discrete domain partition string representation.
 */
std::string DiscreteDomainPartition::asString() const {

    const auto& v = this->asVector();
    std::string s;

    for (size_t i = 0; i < this->v.size(); i++) {
        s += std::to_string(v[i]) + " / ";
    }
    return s.substr(0, s.length() - 3);
}

/**
 *  @return The discrete domain partition vector representation.
 */
std::vector<size_t> DiscreteDomainPartition::asVector() const {
    std::vector<size_t> v;
    for (size_t i = 0; i < this->operator()(0).dimension(); i++) {
        size_t domain = -1;
        for (size_t d = 0; this->dimension(); d++) {
            if (this->operator()(d)(i)) {
                domain = d;
            }
        }

        if (domain == -1) {
            throw std::invalid_argument("DiscreteDomainPartition::asVector(): Found index " + std::to_string(d) + "that does not belong to any domain.");
        } else {
            v.push_back(d);
        }
    }

    return v;
}


}  // namespace GQCP
