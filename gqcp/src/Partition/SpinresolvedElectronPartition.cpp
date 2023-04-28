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


#include "Partition/SpinResolvedElectronPartition.hpp"


namespace GQCP {


/*
 *  MARK: General info
 */

/**
 *  @return The electron partition string representation.
 */
std::string SpinResolvedElectronPartition::asString() const {
    return this->alpha().asString() + std::string(" | ") + this->beta().asString();
}


/**
 *  @return     The number of spin-resolved electrons the partition contains.
 */
size_t SpinResolvedElectronPartition::numberOfElectrons() const {
    return this->alpha().numberOfElectrons() + this->beta().numberOfElectrons();
}


}  // namespace GQCP
