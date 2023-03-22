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


#include "Domain/HubbardDomain.hpp"


namespace GQCP {


/**
 * Calculate the overlap between the Hubbard domain and a spin-unresolved ONV since both can be represented as a bitstring.
 *
 *  @param onv            The spin-unresolved ONV.
 *
 *  @return     The number of overlapping set bits after a bit-by-bit comparison between the Hubbard domain and the spin-unresolved ONV.
 */
size_t HubbardDomain::overlapWithONV(const SpinUnresolvedONV& onv) const {
    boost::dynamic_bitset<> overlap;
    overlap = this->domain & boost::dynamic_bitset<>(this->dimension(), onv.unsignedRepresentation());

    return overlap.count();
}


/**
 * Calculate the overlap between the Hubbard domain and a spin-resolved ONV since the Hubbard domain and each spin-type of the ONV can be represented as a bitstring.
 *
 *  @param onv            The spin-resolved ONV.
 *
 *  @return     The number of overlapping set bits after a bit-by-bit comparison between the Hubbard domain and the spin-types of the spin-resolved ONV.
 */
SpinResolved<size_t> HubbardDomain::overlapWithONV(const SpinResolvedONV& onv) const {
    return SpinResolved<size_t>(this->overlapWithONV(onv.onv(Spin::alpha)), this->overlapWithONV(onv.onv(Spin::beta)));
}


}  // namespace GQCP
