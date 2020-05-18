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
 *  Create an implicit orbital space with the given dimensions.
 * 
 *  @param counts               a map that links an occupation type (k_occupied, k_active, k_virtual) with the number of orbitals that are to be found in that orbital space
 * 
 *  @note An 'implicit' orbital space is one where all indices are sorted by increasing value, and the occupied indices are lower than the active indices, which are in turn lower than the virtual indices.
 */
OrbitalSpace OrbitalSpace::Implicit(const std::map<OccupationType, size_t>& counts) {

    size_t start_index = 0;

    // Create the occupied indices, if any.
    std::vector<size_t> occupied_indices {};
    auto it = counts.find(OccupationType::k_occupied);
    if (it != counts.end()) {
        const auto number_of_occupied_orbitals = it->second;
        occupied_indices = std::vector<size_t>(number_of_occupied_orbitals);
        std::iota(occupied_indices.begin(), occupied_indices.end(), start_index);
        start_index += number_of_occupied_orbitals;
    }

    // Create the active indices, if any.
    std::vector<size_t> active_indices {};
    it = counts.find(OccupationType::k_active);
    if (it != counts.end()) {
        const auto number_of_active_orbitals = it->second;
        active_indices = std::vector<size_t>(number_of_active_orbitals);
        std::iota(active_indices.begin(), active_indices.end(), start_index);
        start_index += number_of_active_orbitals;
    }

    // Create the virtual indices, if any.
    std::vector<size_t> virtual_indices {};
    it = counts.find(OccupationType::k_virtual);
    if (it != counts.end()) {
        const auto number_of_virtual_orbitals = it->second;
        virtual_indices = std::vector<size_t>(number_of_virtual_orbitals);
        std::iota(virtual_indices.begin(), virtual_indices.end(), start_index);
    }


    return OrbitalSpace(occupied_indices, active_indices, virtual_indices);
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
    for (const auto& i : this->indices(OccupationType::k_occupied)) {
        description_string += (boost::format("%s ") % i).str();
    }

    description_string += "\n\tThe active orbital indices are: ";
    for (const auto& m : this->indices(OccupationType::k_active)) {
        description_string += (boost::format("%s ") % m).str();
    }

    description_string += "\n\tThe virtual orbital indices are: ";
    for (const auto& a : this->indices(OccupationType::k_virtual)) {
        description_string += (boost::format("%s ") % a).str();
    }

    return description_string;
}


/**
 *  @param type             the occupation type that the indices should belong to
 * 
 *  @return the indices that belong to the given occupation type
 */
const std::vector<size_t>& OrbitalSpace::indices(const OccupationType type) const {

    switch (type) {
    case OccupationType::k_occupied: {
        return this->occupied_indices;
        break;
    }

    case OccupationType::k_active: {
        return this->active_indices;
        break;
    }

    case OccupationType::k_virtual: {
        return this->virtual_indices;
        break;
    }
    }
}


/**
 *  @param type             an occupation type (k_occupied, k_active, k_virtual)
 *  @param p                an orbital index
 * 
 *  @return if the orbital at the given index is in the given orbital space
 */
bool OrbitalSpace::isIndex(const OccupationType type, const size_t p) const {

    const auto& indices = this->indices(type);
    return std::find(indices.begin(), indices.end(), p) != indices.end();
}


}  // namespace GQCP
