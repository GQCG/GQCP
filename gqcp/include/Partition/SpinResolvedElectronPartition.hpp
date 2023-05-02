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


#include "Partition/SpinUnresolvedElectronPartition.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 * A spin-resolved partition of alpha and beta electron numbers over e.g. domains.
 */
class SpinResolvedElectronPartition:
    public SpinResolved<SpinUnresolvedElectronPartition> {
public:
    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolved<SpinUnresolvedElectronPartition>::Of;

    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolved`'s constructors.
    using SpinResolved<SpinUnresolvedElectronPartition>::SpinResolved;

    /*
     *  MARK: General info
     */

    /**
     *  @return The spin-resolved electron partition string representation.
     */
    std::string asString() const;

    /**
     * @param i     The index in the spin-resolved electron partition.
     *
     *
     * @return     The number of alpha and beta electrons the partition contains at index `i`.
     */
    SpinResolved<size_t> numberOfElectrons(size_t i) const { return SpinResolvedBase<size_t, size_t>(this->alpha().numberOfElectrons(), this->beta().numberOfElectrons()); }

    /**
     *  @return     The number of electrons the partition contains.
     */
    size_t numberOfElectrons() const;
};


}  // namespace GQCP
