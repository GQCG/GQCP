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


#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/**
 *  A type that can be used to build (spin-unresolved) ONV paths in a graphical representation of a CI addressing scheme.
 * 
 *  When modifying paths, instances of this type accordingly update the path's corresponding address and total sign.
 */
class ONVPath {
private:
    // The ONV basis in which the ONV (path) lives. A (const) reference in order to avoid unnecessary copying.
    const SpinUnresolvedONVBasis& onv_basis;

    // The original ONV that encapsulates the indices of the orbitals that are occupied. A (const) reference in order to avoid unnecessary copying.
    const SpinUnresolvedONV& onv;

    // The address of the current path.
    size_t m_address;

    // The orbital index that should be checked next for a possible creation. Since we're always constructing paths from the top-left to the bottom-right, this index will always be larger than the index q on which we previously annihilated. After creation, the path then corresponds to E_{pq} |onv>, with |onv> the initial ONV.
    size_t p;

    // The total phase factor/sign associated to the original path's modification.
    int m_sign;


public:
    // CONSTRUCTORS

    /**
     *  Initialize an ONV path corresponding to the given ONV that lives in the given ONV basis with the given address.
     * 
     *  @param onv_basis            the ONV basis in which the ONV (path) lives
     *  @param onv                  the ONV whose path should be represented
     *  @param address              the address of the ONV inside the given ONV basis
     */
    ONVPath(const SpinUnresolvedONVBasis& onv_basis, const SpinUnresolvedONV& onv, const size_t address);

    /**
     *  Initialize an ONV path corresponding to the given ONV that lives in the given ONV basis.
     * 
     *  @param onv_basis            the ONV basis in which the ONV (path) lives
     *  @param onv                  the ONV whose path should be represented
     */
    ONVPath(const SpinUnresolvedONVBasis& onv_basis, const SpinUnresolvedONV& onv);


    // PUBLIC METHODS

    /**
     *  @return the address of the current path
     */
    size_t address() const { return this->m_address; }

    /**
     *  Annihilate the diagonal arc that starts at the coordinate (q,n).
     * 
     *  @param q        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index q
     */
    void annihilate(const size_t q, const size_t n);

    /**
     *  Create the diagonal arc that starts at the coordinate (p, n).
     * 
     *  @param p        the index of the orbital that should be created
     *  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
     */
    void create(const size_t p, const size_t n);

    /**
     *  Translate the diagonal arc that starts at the coordinate (p, n) to the left.
     * 
     *  @param p        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index p
     */
    void leftTranslate(const size_t p, const size_t n);

    /**
     *  @return The orbital index "p" that should be checked next for a possible creation. Since we're always constructing paths from the top-left to the bottom-right, this index will always be larger than the index "q" on which we previously annihilated. After creation, the path then corresponds to E_{pq} |onv>, with |onv> the initial ONV.
     */
    size_t nextCreationIndex() const { return this->p; }

    /**
     *  @return the total phase factor/sign associated to the original path's modification
     */
    int sign() const { return this->m_sign; }
};


}  // namespace GQCP
