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


//#include "ONVBasis/SpinUnresolvedONV.hpp"
//#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/**
 *  A type that can be used to build (spin-unresolved) ONV paths in a graphical representation of a CI addressing scheme.
 * 
 *  When modifying paths, instances of this type accordingly update the path's corresponding address and total sign.
 */
template <typename ONVBasis>
class ONVPath {
public:
    using ONV = typename ONVBasis::ONV;

private:
    // The ONV basis in which the ONV (path) lives. A (const) reference in order to avoid unnecessary copying.
    const ONVBasis& onv_basis;

    // The original ONV that encapsulates the indices of the orbitals that are occupied. A (const) reference in order to avoid unnecessary copying.
    const ONV& onv;

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
    ONVPath(const ONVBasis& onv_basis, const ONV& onv, const size_t address) :
        onv_basis {onv_basis},
        onv {onv},
        m_address {address},
        orbital_index {0},
        electron_index {0},
        m_sign {+1} {};

    /**
     *  Initialize an ONV path corresponding to the given ONV that lives in the given ONV basis.
     * 
     *  @param onv_basis            the ONV basis in which the ONV (path) lives
     *  @param onv                  the ONV whose path should be represented
     */
    ONVPath(const ONVBasis& onv_basis, const ONV& onv) :
        ONVPath(onv_basis, onv, onv_basis.addressOf(onv)) {};


    // PUBLIC METHODS

    /**
     *  @return the address of the current path
     */
    size_t address() const { return this->m_address; }

    /**
     *  @param p        the index of the orbital that should be created
     *  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
     * 
     *  @return the address of the current path after creation on (p, n).
     */
    size_t addressAfterCreation(const size_t p, const size_t n) {
        return this->m_address + this->onv_basis.arcWeight(p, n);
    }

    /**
     * Create a vertical arc on the path's current state. 
     * 
     * @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void addVertical() { this->orbital_index++; }

    /**
     *  According to this path's current state, annihilate the next diagonal arc.
     *
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void annihilate() {

        // During annihilation, we're removing an arc, so we'll have to update the current address by removing the corresponding arc weight.
        this->m_address -= this->onv_basis.arcWeight(this->orbital_index, this->electron_index);

        // Annihilating the diagonal arc at the current vertex (p,n) means that we will only consider paths where orbital p is unoccupied. Therefore, the sub-path up until vertex (p+1,n) is considered finished.
        this->orbital_index++;

        // Note that we're not calling `this->annihilate(p,n)` in order to avoid a null-operation `this->electron_index = this->electron_index`.
    }

    /**
     *  Annihilate the diagonal arc that starts at the coordinate (q,n).
     * 
     *  @param q        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index q
     * 
     *  @note Call this method only when the given orbital index 'q' is occupied! This method does not perform any validation checks.
     */
    void annihilate(const size_t q, const size_t n) {

        // During annihilation, we're removing an arc, so we'll have to update the current address by removing the corresponding arc weight.
        this->m_address -= this->onv_basis.arcWeight(q, n);

        // Annihilating the diagonal arc at (q,n) means that we will only consider paths where orbital q is unoccupied. Therefore, the sub-path up until vertex (q+1,n) is considered finished.
        this->orbital_index = q + 1;
        this->electron_index = n;
    }

    /**
     *  According to this path's current state, create the next diagonal arc.
     * 
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void create() {

        this->create(this->orbitalIndex(), this->electronIndex());

        // Note that, contrary to `this->annihilate()`, we may call `this->create(p,n)` because the number of operations to be done is the same.
    }

    /**
     *  Create the diagonal arc that starts at the coordinate (p, n).
     * 
     *  @param p        the index of the orbital that should be created
     *  @param n        the number of electrons in the ONV/path up to the orbital index q, prior to the creation
     */
    void create(const size_t p, const size_t n) {

        // During creation, we're adding an arc, so we'll have to update the current address by adding the corresponding arc weight.
        this->m_address += this->onv_basis.arcWeight(p, n);

        // Creating the diagonal arc at (p,n) means that we will only consider paths where orbital p is occupied. Therefore, the sub-path up until vertex (p+1,n+1) is considered finished.
        this->orbital_index = p + 1;
        this->electron_index = n + 1;
    }


    /**
     *  @return The electron index 'n' that, together with the orbital index 'p' signifies the vertex (p,n) up until which the ONV path construction is finished.
     *
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    size_t electronIndex() const { return this->electron_index; }

    /**
     *  @return If the path's construction is considered finished.
     */
    bool isFinished() const {
        return this->electron_index >= this->onv_basis.numberOfElectrons() || this->orbital_index >= this->onv_basis.numberOfOrbitals();
    }

    /**
     *  @return If the orbital index 'p' is allowed. If p exceeds a certain value, there are not enough electrons for the given set of orbitals.
     */
    bool isOrbitalIndexValid() const {
        return this->orbitalIndex() <= (this->onv_basis.numberOfOrbitals() - this->onv_basis.numberOfElectrons() + this->electronIndex());
    }

    /**
     *  Translate the diagonal arc that starts at the coordinate (p,n) to the left, indicating that the current path is 'open' at the vertex (p,n-1) and that the orbital 'p' should be occupied in subsequent path manipulations.
     * 
     *  @param p        the index of the orbital that should be annihilated
     *  @param n        the number of electrons in the ONV/path up to the orbital index p
     * 
     *  @note Call this method only if you're sure that the path is 'open' at the vertex (p, n-1)! This method does not perform any validation checks.
     */
    void leftTranslate(const size_t p, const size_t n) {

        // In order to keep the operation count as small as possible, we're not using `this->create()` or `this->annihilate()`.

        // Translating a diagonal arc can be rewritten as a removal of the arc weight (p,n), followed by the addition of the arc weight (p,n-1).
        this->m_address += this->onv_basis.arcWeight(p, n - 1) - this->onv_basis.arcWeight(p, n);

        // Since a left-translation describes the process of 'encountering an electron/occupied orbital', the sign factor should be updated according to the fermionic anticommutation rules.
        this->m_sign *= -1;

        // In subsequent path manipulations, the orbital `p` should be occupied (while no 'net creation' has occurred), so the sub path's construction is considered fixed up until the indices (p+1,n).
        this->orbital_index = p + 1;
        this->electron_index = n;
    }

    /**
     *  According to this path's current state, translate diagonal arcs to the left until an unoccupied orbital (vertical arc) is found.
     * 
     *  @note The current state of this path can be retrieved by inspecting the point (p,n), where `p = this->orbitalIndex()` and `n = this->electronIndex()`, which signifies the vertex up until which the path construction is complete.
     */
    void leftTranslateUntilVertical() {

        while (!this->isFinished() && this->onv.isOccupied(this->orbital_index)) {

            // Translate the diagonal arc starting at (p,n+1) one position to the left. This function keeps tabs on the creation index, electron index and sign.
            this->leftTranslate(this->orbital_index, this->electron_index + 1);
        }
    }

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
