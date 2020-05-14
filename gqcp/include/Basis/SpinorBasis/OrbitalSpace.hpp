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


#include <cstdlib>
#include <string>
#include <vector>


namespace GQCP {


/**
 *  A class that encapsulates occupied, active and virtual orbital indices.
 * 
 *  @note The union of these three sets of indices is supposed to be the full set of orbital indices.
 *  @note This class is intended to be used in conjunction with restricted spin-orbital bases (RSpinorBasis) or general spinor bases (GSpinorBasis).
 */
class OrbitalSpace {
private:
    std::vector<size_t> occupied_indices;  // the indices of the orbitals that are considered occupied by the electrons
    std::vector<size_t> active_indices;    // the indices of the orbitals that are considered 'active', i.e. those spinor indices that are occupied in some set of configurations but not in others
    std::vector<size_t> virtual_indices;   // the indices of the orbitals that are considered virtual (i.e. unoccupied) by the electrons

    std::vector<size_t> all_indices;


public:
    // CONSTRUCTORS

    /**
     *  Construct an orbital space from occupied, active and virtual indices.
     * 
     *  @param occupied_indices                 the indices of the orbitals that are considered occupied by the electrons
     *  @param active_indices                   the indices of the orbitals that are considered 'active', i.e. those spinor indices that are occupied in some set of configurations but not in others
     *  @param virtual_indices                  the indices of the orbitals that are considered virtual (i.e. unoccupied) by the electrons
     */
    OrbitalSpace(const std::vector<size_t>& occupied_indices, const std::vector<size_t>& active_indices, const std::vector<size_t>& virtual_indices);

    /**
     *  Construct an orbital space only from occupied and virtual indices.
     * 
     *  @param occupied_indices                 the indices of the orbitals that are considered occupied by the electrons
     *  @param virtual_indices                  the indices of the orbitals that are considered virtual (i.e. unoccupied) by the electrons
     * 
     *  @note This constructor assumes that there are no active orbital indices
     */
    OrbitalSpace(const std::vector<size_t>& occupied_indices, const std::vector<size_t>& virtual_indices);


    // NAMED CONSTRUCTORS

    /**
     *  Create an orbital space with only occupied indices [0, N[.
     * 
     *  @param N                the number of occupied indices
     * 
     *  @return an orbital space with only occupied indices [0, N[.
     */
    static OrbitalSpace Occupied(const size_t N);

    /**
     *  Create an orbital space that is separated between occupied and virtual orbitals.
     * 
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     * 
     *  @return an orbital space that is separated between occupied and virtual orbitals.
     */
    static OrbitalSpace OccupiedVirtual(const size_t N, const size_t M);


    // PUBLIC METHODS

    /**
     *  @return all the indices of the spinors
     */
    const std::vector<size_t>& allIndices() const { return this->all_indices; }

    /**
     *  @return the active orbital space, i.e. those spinor indices that are occupied in some set of configurations but not in others
     */
    const std::vector<size_t>& activeIndices() const { return this->active_indices; }

    /**
     *  @return a textual description of this orbital space
     */
    std::string description() const;

    /**
     *  @param p            an orbital index
     * 
     *  @return if the orbital at the given index is in the active orbital space
     */
    bool isIndexActive(const size_t p) const;

    /**
     *  @param p            an orbital index
     * 
     *  @return if the orbital at the given index is in the occupied orbital space
     */
    bool isIndexOccupied(const size_t p) const;

    /**
     *  @param p            an orbital index
     * 
     *  @return if the orbital at the given index is in the virtual orbital space
     */
    bool isIndexVirtual(const size_t p) const;

    /**
     *  @return the total number of orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space
     */
    size_t numberOfOrbitals() const { return this->allIndices().size(); }

    /**
     *  @return the number of active orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space
     */
    size_t numberOfActiveOrbitals() const { return this->activeIndices().size(); }

    /**
     *  @return the number of occupied orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space
     */
    size_t numberOfOccupiedOrbitals() const { return this->occupiedIndices().size(); }

    /**
     *  @return the number of virtual orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space
     */
    size_t numberOfVirtualOrbitals() const { return this->virtualIndices().size(); }

    /**
     *  @return the occupied orbital space, i.e. the indices of the orbitals that are considered occupied by the electrons
     */
    const std::vector<size_t>& occupiedIndices() const { return this->occupied_indices; }

    /**
     *  @return the virtual orbital space, i.e. the indices of the orbitals that are considered virtual by the electrons
     */
    const std::vector<size_t>& virtualIndices() const { return this->virtual_indices; }
};


}  // namespace GQCP
