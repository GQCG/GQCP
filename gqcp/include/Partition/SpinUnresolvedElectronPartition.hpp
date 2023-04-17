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

#include <numeric>


namespace GQCP {


/**
 * A partition of electron numbers over e.g. domains.
 */
class SpinUnresolvedElectronPartition:
    public SimplePartition<SpinUnresolvedElectronPartition> {

public:
    /*
     *  MARK: General info
     */

    /**
     *  @return The electron partition string representation.
     */
    std::string asString() const;

    /**
     * @param i     The index in the electron partition.
     *
     *
     * @return     The number of electrons the partition contains at index `i`.
     */
    size_t numberOfElectrons(size_t i) const { return this->operator()(i); }

    /**
     *  @return     The number of electrons the partition contains.
     */
    size_t numberOfElectrons() const;
};


}  // namespace GQCP
