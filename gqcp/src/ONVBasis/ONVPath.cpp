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

#include "ONVBasis/ONVPath.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  Initialize an ONV path corresponding to the given ONV that lives in the given ONV basis with the given address.
 * 
 *  @param onv_basis            the ONV basis in which the ONV (path) lives
 *  @param onv                  the ONV whose path should be represented
 *  @param address              the address of the ONV inside the given ONV basis
 */
ONVPath::ONVPath(const SpinUnresolvedONVBasis& onv_basis, const SpinUnresolvedONV& onv, const size_t address) :
    onv_basis {onv_basis},
    onv {onv},
    m_address {address},
    p {0},
    m_sign {+1} {}


/**
 *  Initialize an ONV path corresponding to the given ONV that lives in the given ONV basis.
 * 
 *  @param onv_basis            the ONV basis in which the ONV (path) lives
 *  @param onv                  the ONV whose path should be represented
 */
ONVPath::ONVPath(const SpinUnresolvedONVBasis& onv_basis, const SpinUnresolvedONV& onv) :
    ONVPath(onv_basis, onv, onv_basis.addressOf(onv)) {}


/*
 *  PUBLIC METHODS
 */

/**
 *  Annihilate the diagonal arc that starts at the coordinate (q,n).
 * 
 *  @param q        the index of the orbital that should be annihilated
 *  @param n        the number of electrons in the ONV/path up to the orbital index q
 */
void ONVPath::annihilate(const size_t q, const size_t n) {

    // During annihilation, we're removing an arc, so we'll have to update the current address by removing the corresponding arc weight.
    this->m_address -= this->onv_basis.arcWeight(q, n);

    // Update the next possible creation index. Since we're always constructing paths from the top-left to the bottom-right, we're only considering creation indices p > q.
    this->p = q + 1;
}


/**
 *  Create the diagonal arc that starts at the coordinate (p, n).
 * 
 *  @param p        the index of the orbital that should be created
 *  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
 */
void ONVPath::create(const size_t p, const size_t n) {

    // During creation, we're adding an arc, so we'll have to update the current address by adding the corresponding arc weight.
    this->m_address += this->onv_basis.arcWeight(p, n);
}


/**
 *  Translate the diagonal arc that starts at the coordinate (p, n) to the left.
 * 
 *  @param p        the index of the orbital that should be annihilated
 *  @param n        the number of electrons in the ONV/path up to the orbital index p
 */
void ONVPath::leftTranslate(const size_t p, const size_t n) {

    // Translating a diagonal arc can be rewritten as an annihilation, followed by a creation.
    this->annihilate(p, n);
    this->create(p, n - 1);

    // Since a left-translation describes the process of 'encountering an electron/occupied orbital', the sign factor should be updated according to the fermionic anticommutation rules.
    this->m_sign *= -1;
}

/**
 * Close the open path by shifting diagonal arcs to the left. Stop when an unoccupied orbital (vertical arc) is found.
 * 
 * @param p         index of the orbital from where we start closing the path
 * @param n         the number of electrons in the ONV/path up to the orbital index p
 */
void ONVPath::shiftUntilNextUnoccupiedOrbital(size_t n) {
    
    // If the orbital index is not the same as the occupation index of the next electron, we have encountered an unoccupied orbital/vertical arc.
    while (n < this->onv_basis.N && this->p == this->onv.occupationIndexOf(n+1)) {

        // Translate the diagonal arc starting at (p,n+1) one position to the left. This function keeps tabs on the creation index p and sign.
        this->leftTranslate(this->p, n+1);
        
        // Shift to next electron to be checked in the loop.
        n++;
    }
}


}  // namespace GQCP
