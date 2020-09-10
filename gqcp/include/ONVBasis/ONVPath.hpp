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

    // The orbital index 'p' that, together with the electron index 'n' signifies the vertex (p,n) up until which the ONV path construction is finished.
    size_t orbital_index;

    // The electron index 'n' that, together with the orbital index 'p' signifies the vertex (p,n) up until which the ONV path construction is finished.
    size_t electron_index;

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
     *  According to this path's current state, annihilate the next diagonal arc.
     *
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void annihilate();

    /**
     *  Annihilate the diagonal arc that starts at the coordinate (q,n).
     * 
     *  @param q        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index q
     * 
     *  @note Call this method only when the given orbital index 'q' is occupied! This method does not perform any validation checks.
     */
    void annihilate(const size_t q, const size_t n);

    /**
     *  According to this path's current state, create the next diagonal arc.
     * 
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void create();

    /**
     *  Create the diagonal arc that starts at the coordinate (p, n).
     * 
     *  @param p        the index of the orbital that should be created
     *  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
     */
    void create(const size_t p, const size_t n);

    /**
     *  @return The electron index 'n' that, together with the orbital index 'p' signifies the vertex (p,n) up until which the ONV path construction is finished.
     *
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    size_t electronIndex() const { return this->electron_index; }

    /**
     *  @return If the path's construction is considered finished.
     */
    bool isFinished() const;

    /**
     *  Translate the diagonal arc that starts at the coordinate (p,n) to the left, indicating that the current path is 'open' at the vertex (p,n-1) and that the orbital 'p' should be occupied in subsequent path manipulations.
     * 
     *  @param p        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index p
     * 
     *  @note Call this method only if you're sure that the path is 'open' at the vertex (p, n-1)! This method does not perform any validation checks.
     */
    void leftTranslate(const size_t p, const size_t n);

    /**
     *  According to this path's current state, translate diagonal arcs to the left until an unoccupied orbital (vertical arc) is found.
     * 
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void leftTranslateUntilVertical();

    /**
     *  @return The orbital index 'p' that, together with the electron index 'n' signifies the vertex (p,n) up until which the ONV path construction is finished.
     *
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    size_t orbitalIndex() const { return this->orbital_index; }

    /**
     *  @return the total phase factor/sign associated to the original path's modification
     */
    int sign() const { return this->m_sign; }
};


}  // namespace GQCP
