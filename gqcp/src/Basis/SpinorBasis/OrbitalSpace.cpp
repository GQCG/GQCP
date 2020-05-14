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

#include "Basis/SpinorBasis/OrbitalSpace.hpp"

#include <boost/format.hpp>

#include <numeric>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Construct an orbital space from occupied, active and virtual indices.
 * 
 *  @param occupied_indices                 the indices of the orbitals that are considered occupied by the electrons
 *  @param active_indices                   the indices of the orbitals that are considered 'active', i.e. those spinor indices that are occupied in some set of configurations but not in others
 *  @param virtual_indices                  the indices of the orbitals that are considered virtual (i.e. unoccupied) by the electrons
 */
OrbitalSpace::OrbitalSpace(const std::vector<size_t>& occupied_indices, const std::vector<size_t>& active_indices, const std::vector<size_t>& virtual_indices) :
    occupied_indices {occupied_indices},
    active_indices {active_indices},
    virtual_indices {virtual_indices} {

    // Generate the list of all indices by taking all elements from the occupied, active and virtual indices.
    this->all_indices = std::vector<size_t>();

    this->all_indices.insert(all_indices.end(), occupied_indices.begin(), occupied_indices.end());
    this->all_indices.insert(all_indices.end(), active_indices.begin(), active_indices.end());
    this->all_indices.insert(all_indices.end(), virtual_indices.begin(), virtual_indices.end());
}


/**
 *  Construct an orbital space from occupied and virtual indices.
 * 
 *  @param occupied_indices                 the indices of the orbitals that are considered occupied by the electrons
 *  @param virtual_indices                  the indices of the orbitals that are considered virtual (i.e. unoccupied) by the electrons
 */
OrbitalSpace::OrbitalSpace(const std::vector<size_t>& occupied_indices, const std::vector<size_t>& virtual_indices) :
    OrbitalSpace(occupied_indices, {}, virtual_indices) {}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Create an orbital space with only occupied indices [0, N[.
 * 
 *  @param N                the number of occupied indices
 * 
 *  @return an orbital space with only occupied indices [0, N[.
 */
OrbitalSpace OrbitalSpace::Occupied(const size_t N) {

    std::vector<size_t> occupied_indices(N);  // zero-initialized with N entries

    std::iota(occupied_indices.begin(), occupied_indices.end(), 0);  // fill with 0 to N-1
    return OrbitalSpace(occupied_indices, {}, {});
}


/**
 *  Create an orbital space that is separated between occupied and virtual orbitals.
 * 
 *  @param N                the number of occupied orbitals
 *  @param M                the total number of orbitals
 * 
 *  @return an orbital space that is separated between occupied and virtual orbitals.
 */
OrbitalSpace OrbitalSpace::OccupiedVirtual(const size_t N, const size_t M) {

    std::vector<size_t> occupied_indices(N);     // zero-initialized with N entries
    std::vector<size_t> virtual_indices(M - N);  // zero-initialized with (M-N) entries

    std::iota(occupied_indices.begin(), occupied_indices.end(), 0);  // fill with 0 to N-1
    std::iota(virtual_indices.begin(), virtual_indices.end(), N);    // fill with N to K-1

    return OrbitalSpace(occupied_indices, virtual_indices);
}


/*
 *  PUBLIC METHODS
 */


/**
 *  @return a textual description of this orbital space
 */
std::string OrbitalSpace::description() const {

    std::string description_string = (boost::format("An orbital space with %s orbitals.\n") % this->numberOfOrbitals()).str();

    description_string += "\n\tThe occupied orbital indices are: ";
    for (const auto& i : this->occupiedIndices()) {
        description_string += (boost::format("%s ") % i).str();
    }

    description_string += "\n\tThe active orbital indices are: ";
    for (const auto& m : this->activeIndices()) {
        description_string += (boost::format("%s ") % m).str();
    }

    description_string += "\n\tThe virtual orbital indices are: ";
    for (const auto& a : this->virtualIndices()) {
        description_string += (boost::format("%s ") % a).str();
    }

    return description_string;
}


/**
 *  @param p            an orbital index
 * 
 *  @return if the orbital at the given index is in the active orbital space
 */
bool OrbitalSpace::isIndexActive(const size_t p) const {

    return std::find(this->activeIndices().begin(), this->activeIndices().end(), p) != this->activeIndices().end();
}


/**
 *  @param p            an orbital index
 * 
 *  @return if the orbital at the given index is in the occupied orbital space
 */
bool OrbitalSpace::isIndexOccupied(const size_t p) const {

    return std::find(this->occupiedIndices().begin(), this->occupiedIndices().end(), p) != this->occupiedIndices().end();
}


/**
 *  @param p            an orbital index
 * 
 *  @return if the orbital at the given index is in the virtual orbital space
 */
bool OrbitalSpace::isIndexVirtual(const size_t p) const {

    return std::find(this->virtualIndices().begin(), this->virtualIndices().end(), p) != this->virtualIndices().end();
}


}  // namespace GQCP
