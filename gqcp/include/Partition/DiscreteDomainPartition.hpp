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


/**
 * A partition (i.e., collection) of discrete domains.
 */
class DiscreteDomainPartition:
    public DomainPartition<DiscreteDomain> {

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a discrete domain partition from a vector of discrete domains.
     *
     *  @param domains          The vector of discrete domains.
     */
    DiscreteDomainPartition(const std::vector<DiscreteDomain>& domains);

    /**
     *  Create a discrete domain partition from a vector of unsigned representations.
     *
     *  @param domains          The vector of unsigned representations.
     *  @param M                The dimension of each domain.
     */
    DiscreteDomainPartition(const std::vector<size_t>& unsigned_representations, size_t M);

    /*
     *  MARK: General info
     */

    /**
     *  @return The discrete domain partition string representation.
     */
    std::string asString() const;

    /**
     *  @return The discrete domain partition vector representation.
     */
    std::vector<size_t> asVector() const;
};


}  // namespace GQCP
