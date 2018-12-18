// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef GQCP_FOCKSPACE_HPP
#define GQCP_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"



namespace GQCP {


/**
 *  The full Fock space for a number of orbitals and number of electrons
 *
 *  The ONVs and addresses are linked with a hashing function calculated with an addressing scheme. The implementation of the addressing scheme is from Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 *
 */
class FockSpace: public BaseFockSpace {
private:
    size_t N;  // number of electrons
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    // PRIVATE METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation) const;


public:
    // CONSTRUCTORS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    const Matrixu& get_vertex_weights() const { return this->vertex_weights; }
    size_t get_N() const { return this->N; }
    FockSpaceType get_type() const override { return FockSpaceType::FockSpace; }


    // ITERATOR
    class Iterator {
    private:
        const FockSpace* fock_space;  // non-owning pointer: the Fock space is never in a stack level deeper than its own iterator
        size_t address;  // the address of the current onv
        ONV onv;  // the current onv


        // CONSTRUCTORS
        /**
         *  Place an iterator in a possibly invalid state (the ONV does not correspond to the address)
         *
         *  @note This constructor is only implemented to provide a proper FockSpace.end() iterator with an address that is higher than the 'last' ONV in the Fock space
         *
         *  @param fock_space       the Fock space that should be iterated over
         *  @param address          the address of the current ONV
         *  @param onv              the current ONV
         */
        Iterator(const FockSpace& fock_space, size_t address, const ONV& onv);
        friend class FockSpace;  // declare the whole class friend because we don't want to forward declare Iterator and then make the FockSpace::end() function a friend


    public:
        // CONSTRUCTORS
        /**
         *  @param fock_space       the Fock space that should be iterated over
         *  @param address          the address of the current ONV
         */
        Iterator(const FockSpace& fock_space, size_t address);

        /**
         *  Constructor that starts at the ONV with address 0
         *
         *  @param fock_space       the Fock space that should be iterated over
         */
        Iterator(const FockSpace& fock_space);


        // OPERATORS
        /**
         *  Move the iterator forward (pre-increment)
         *
         *  @return a reference to the updated iterator
         */
        Iterator& operator++();

        /**
         *  @param other            the other iterator
         *
         *  @return if this iterator is the same as the other
         */
        bool operator==(const Iterator& other) const;

        /**
         *  @param other            the other iterator
         *
         *  @return if this iterator is not the same as the other
         */
        bool operator!=(const Iterator& other) const;

        /**
         *  @return the current ONV
         */
        const ONV& operator*() const;


        // PUBLIC METHODS
        /**
         *  @return the current ONV
         */
        const ONV& currentONV() const;

        /**
         *  @return the current address
         */
        size_t currentAddress() const;
    };


    // STATIC PUBLIC METHODS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     *
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);


    // PUBLIC METHODS
    /**
     *  @param address      the address (i.e. the ordening number) of the ONV
     *
     *  @return the ONV with the corresponding address
     */
    ONV makeONV(size_t address) const;

    /**
     *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
     *
     *  @param onv      the current ONV
     */
    void setNextONV(ONV& onv) const;

    /**
     *  @param onv      the ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(const ONV& onv) const;
  
    /**
     *  Transform an ONV to one corresponding to the given address
     *
     *  @param onv          the ONV
     *  @param address      the address to which the ONV will be set
     */
    void transformONV(ONV& onv, size_t address) const;

    /**
     *  Find the next unoccupied orbital in a given ONV,
     *  update the electron count, orbital index,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the amount of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     */
    template<int T>
    void shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& address, size_t& q, size_t& e) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.get_occupation_index(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->get_vertex_weights(q, e + 1 - T) - this->get_vertex_weights(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
        }
    }

    /**
     *  Find the next unoccupied orbital in a given ONV,
     *  update the electron count, orbital index, sign,
     *  and update the address by calculating a shift
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @tparam T        the amount of previously annihilated electrons
     *
     *  @param address   the address which is updated
     *  @param onv       the ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param sign      the sign which is flipped for each iteration
     */
    template<int T>
    void shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& address, size_t& q, size_t& e, int& sign) const {

        // Test whether the current orbital index is occupied
        while (e < this->N && q == onv.get_occupation_index(e)) {

            // Take the difference of vertex weights for the encountered electron weights to that of a vertex weight path with "a" fewer electrons
            // +1 is added to the electron index, because of how the addressing scheme is arrayed.
            address += this->get_vertex_weights(q, e + 1 - T) - this->get_vertex_weights(q, e + 1);

            // move to the next electron and orbital
            e++;
            q++;
            sign *= -1;
        }
    }

    /**
     *  @return an iterator pointing at the first ONV in the Fock space
     */
    Iterator begin();

    /**
     *  @return an iterator pointing at the last ONV in the Fock space
     */
    Iterator end();
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
