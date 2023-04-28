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


#include "Partition/SimplePartition.hpp"


namespace GQCP {


/**
 * A partition (i.e., collection) of domains.
 */
template <typename DomainType>
class DomainPartition:
    public SimplePartition<DomainPartition<DomainType>> {

public:
    /*
     *  MARK: General info
     */

    /**
     *  @return The domain partition string representation.
     */
    std::string asString() const {
        std::string s;

        for (size_t i = 0; i < this->dimension(); i++) {
            s += this->operator()(i).asString() + " / ";
        }
        return s.substr(0, s.length() - 3);
    }
};


}  // namespace GQCP
