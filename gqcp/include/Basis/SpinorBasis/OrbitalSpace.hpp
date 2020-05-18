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


#include "Basis/SpinorBasis/OccupationType.hpp"
#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"

#include <cstdlib>
#include <map>
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
     *  Create an implicit orbital space with the given dimensions.
     * 
     *  @param counts               a map that links an occupation type (k_occupied, k_active, k_virtual) with the number of orbitals that are to be found in that orbital space
     * 
     *  @note An 'implicit' orbital space is one where all indices are sorted by increasing value, and the occupied indices are lower than the active indices, which are in turn lower than the virtual indices.
     */
    static OrbitalSpace Implicit(const std::map<OccupationType, size_t>& counts);


    // PUBLIC METHODS

    /**
     *  @return a textual description of this orbital space
     */
    std::string description() const;

    /**
     *  @return all the indices of the spinors
     */
    const std::vector<size_t>& indices() const { return this->all_indices; }

    /**
     *  @param type             the occupation type that the indices should belong to
     * 
     *  @return the indices that belong to the given occupation type
     */
    const std::vector<size_t>& indices(const OccupationType type) const;

    /**
     *  @param type             an occupation type (k_occupied, k_active, k_virtual)
     *  @param p                an orbital index
     * 
     *  @return if the orbital at the given index is in the given orbital space
     */
    bool isIndex(const OccupationType type, const size_t p) const;

    /**
     *  @return the total number of orbitals (i.e. spatial orbitals or spinors, depending on the context) in this orbital space
     */
    size_t numberOfOrbitals() const { return this->indices().size(); }

    /**
     *  @param type             an occupation type (k_occupied, k_active, k_virtual)
     * 
     *  @return the number of orbitals (i.e. spatial orbitals or spinors, depending on the context) that belong to the given occupation type
     */
    size_t numberOfOrbitals(const OccupationType type) const { return this->indices(type).size(); }

    /**
     *  Create a zero-initialized mathematical object that can serve as the representation of a object with the given occupation types.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit matrix
     * 
     *  @param row_type                     the spinor occupation type for the rows
     *  @param column_type                  the spinor occupation type for the columns
     * 
     *  @return an implicit matrix slice, according to the given occupation types
     * 
     *  @note For the representation of an occupied-virtual object (for example the T1-coupled-cluster amplitudes t_i^a), the following method can be called
     *      orbital_space.representableObject(OccupationType::occupied, OccupationType::virtual)
     */
    template <typename Scalar>
    ImplicitMatrixSlice<Scalar> representableObject(const OccupationType row_type, const OccupationType column_type) const {

        // Prepare the necessary members for ImplicitMatrixSlice.
        const auto row_indices = this->indices(row_type);
        const auto column_indices = this->indices(column_type);

        const auto rows = row_indices.size();
        const auto columns = column_indices.size();
        const MatrixX<Scalar> M = MatrixX<Scalar>::Zero(rows, columns);

        return ImplicitMatrixSlice<Scalar>(rows, columns, M);
    }


    /**
     *  Create a zero-initializedmathematical object that can serve as the representation of a object with the given occupation types.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit matrix
     * 
     *  @param axis1_type                  the spinor occupation type for the first tensor axis
     *  @param axis2_type                  the spinor occupation type for the second tensor axis
     *  @param axis3_type                  the spinor occupation type for the third tensor axis
     *  @param axis4_type                  the spinor occupation type for the fourth tensor axis
     * 
     *  @return an implicit rank-four tensor slice, according to the given occupation types
     * 
     *  @note For the representation of an occupied-virtual-occupied-virtual object (for example the T2-coupled-cluster amplitudes t_{ij}^{ab}), the following method can be called
     *      orbital_space.representableObject(OccupationType::occupied, OccupationType::virtual, OccupationType::occupied, OccupationType::virtual)
     */
    template <typename Scalar>
    ImplicitRankFourTensorSlice<Scalar> representableObject(const OccupationType axis1_type, const OccupationType axis2_type, const OccupationType axis3_type, const OccupationType axis4_type) const;
};


}  // namespace GQCP
