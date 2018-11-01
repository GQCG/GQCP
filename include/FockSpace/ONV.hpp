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
#ifndef GQCP_ONV_HPP
#define GQCP_ONV_HPP


#include <Eigen/Dense>

#include "common.hpp"



namespace GQCP {


/**
 *      An ONV in quantum chemistry is a string of creation operators acting on top of a vacuum state.
 *      An example for 3 alpha electrons in a Fock space spanned by 4 spatial orbitals is
 *          a_1^\dagger a_2^\dagger a_3^\dagger |vac> = |1,1,1,0>
 *
 *      In this code bitstrings are read from right to left. This means that the
 *      least significant bit relates to the first orbital. Using this notation is how normally bits are read, leading
 *      to more efficient code. As is also usual, the least significant bit has index 0.
 *          The previous example is then represented by the bit string "0111" (7).
 */
class ONV {
private:
    size_t K;  // number of spatial orbitals
    size_t N;  // number of electrons
    size_t unsigned_representation;  // unsigned representation
    VectorXs occupation_indices;  // the occupied orbital electron indexes
                                  // it is a vector of N elements in which occupation_indices[j]
                                  // gives the occupied orbital index for electron j


public:
    // CONSTRUCTORS
    ONV() = default;
    /**
     *  Constructor
     *  @param K a given number of orbitals
     *  @param N a given number of electrons
     *  @param unsigned_representation a representation for the ONV
     */
    ONV(size_t K, size_t N, size_t unsigned_representation);

    /**
     *  Constructor
     *  @param K a given number of orbitals
     *  @param unsigned_representation a representation for the ONV
     */
    ONV(size_t K, size_t unsigned_representation);


    // OPERATORS
    /**
     *  Overloading of operator<< for a GQCP::ONV to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::ONV& onv);

    /**
     *  @return if this->unsigned_representation equals @param other.unsigned_representation
     */
    bool operator==(ONV& other) const;

    /**
     *  @return if this->unsigned_representation does not equal @param other.unsigned_representation
     */
    bool operator!=(ONV& other) const;


    // GETTERS & SETTERS
    /**
     *  @set to a new representation and calls this->updateOccupationIndices()
     */
    void set_representation(size_t unsigned_representation);
    size_t get_unsigned_representation() const { return unsigned_representation; }
    const VectorXs& get_occupation_indices() const { return occupation_indices; }

    /**
     *  @return index of occupied orbital based on the @param electron index
     *
     */
    size_t get_occupied_index(size_t electron_index) const { return occupation_indices(electron_index); }


    // PUBLIC METHODS
    /**
     *  Extracts the positions of the set bits from the this->unsigned_representation
     *  and places them in the this->occupation_indices
     */
    void updateOccupationIndices();

    /**
     *  @return if the @param p-th spatial orbital is occupied, starting from 0
     *  @param p is counted from right to left
     */
    bool isOccupied(size_t p) const;

    /**
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital
     *  Subsequently perform an in-place annihilation on the orbital @param p
     *
     *  !!! IMPORTANT: does not update this->occupation_indices if required call this->updateOccupationIndices !!!
     */
    bool annihilate(size_t p);

    /**
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital
     *  Subsequently perform an in-place annihilation on the orbital @param p
     *
     *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
     *
     *  !!! IMPORTANT: does not update this->occupation_indices if required call this->updateOccupationIndices !!!
     */
    bool annihilate(size_t p, int& sign);

    /**
     * @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital
     * Subsequently perform an in-place creation on the orbital @param p
     *
     *  !!! IMPORTANT: does not update this->occupation_indices if required call this->updateOccupationIndices !!!
     */
    bool create(size_t p);

    /**
     *  @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital
     *  Subsequently perform an in-place creation on the orbital @param p
     *
     *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
     *
     *  !!! IMPORTANT: does not update this->occupation_indices if required call this->updateOccupationIndices !!!
     */
    bool create(size_t p, int& sign);

    /**
     *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on orbital @param p, starting from 0, read from right to left.
     *
     *  Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
     */
    int operatorPhaseFactor(size_t p) const;


    /**
     *  @return the representation of a slice (i.e. a subset) of the spin string between @param index_start (included)
     *  and @param index_end (not included).
     *
     *  Both @param index_start and @param index_end are read from right to left, which means that the slice
     *  is from right to left as well.
     *
     *      Example:
     *          "010011".slice(1, 4) => "01[001]1" -> "001"
     *
     */
    size_t slice(size_t index_start, size_t index_end) const;


    /**
     *  @return the number of different bits between this and @param other, i.e. two times the number of electron excitations
     */
    size_t countNumberOfDifferences(const ONV& other) const;


    /**
     *  @return the positions of the bits (from right to left) that are occupied in this, but unoccupied in @param other
     */
    std::vector<size_t> findDifferentOccupations(const ONV &other) const;


    /**
     *  @return the positions of the bits (from right to left) that are occupied in @this and occupied in @param other
     */
    std::vector<size_t> findMatchingOccupations(const ONV& other) const;


    /**
     * @return std::string containing the ONV representation
     */
    std::string asString() const;
};


}  // namespace GQCP

#endif  // GQCP_ONV_HPP
