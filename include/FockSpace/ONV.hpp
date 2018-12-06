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
#include <tuple>



namespace GQCP {


/**
 *  A class that represents an ONV (occupation number vector)

 *  An ONV in quantum chemistry is a string of creation operators acting on top of a vacuum state.
 *  An example for 3 alpha electrons in a Fock space spanned by 4 spatial orbitals is
 *      a_1^\dagger a_2^\dagger a_3^\dagger |vac> = |1,1,1,0>
 *
 *  In this code bitstrings are read from right to left. This means that the least significant bit relates to the first orbital.
 *  Using this notation is how normally bits are read, leading to more efficient code.
 *  As is also usual, the least significant bit has index 0. The previous example is then represented by the bit string "0111" (7).
 */
class ONV {
private:
    size_t K;  // number of spatial orbitals
    size_t N;  // number of electrons
    size_t unsigned_representation;
    VectorXs occupation_indices;  // the occupied orbital electron indices
                                  // it is a vector of N elements in which occupation_indices[j]
                                  // gives the occupied orbital index for electron j


public:
    // CONSTRUCTORS
    ONV() = default;
    /**
     *  @param K                        the number of orbitals
     *  @param N                        the number of electrons
     *  @param unsigned_representation  the representation for the ONV as an unsigned integer
     */
    ONV(size_t K, size_t N, size_t unsigned_representation);

    /**
     *  Constructs a default ONV without a representation
     *
     *  @param K                        the number of orbitals
     *  @param N                        the number of electrons
     */
    ONV(size_t K, size_t N);


    // OPERATORS
    /**
     *  Overloading of operator<< for a GQCP::ONV to be used with streams
     *
     *  @param os       the output stream which the ONV should be concatenated to
     *  @param onv      the ONV that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::ONV& onv);

    /**
     *  @param other    the other ONV
     *
     *  @return if this ONV is the same as the other ONV
     */
    bool operator==(ONV& other) const;

    /**
     *  @param other    the other ONV
     *
     *  @return if this ONV is not the same as the other ONV
     */
    bool operator!=(ONV& other) const;


    // SETTERS
    /**
     *  @param unsigned_representation      the new representation as an unsigned integer
     *
     *  Set the representation of an ONV to a new representation and call update the occupation indices accordingly
     */
    void set_representation(size_t unsigned_representation);


    // GETTERS
    size_t get_unsigned_representation() const { return unsigned_representation; }
    const VectorXs& get_occupation_indices() const { return occupation_indices; }

    /**
     *  @param electron_index       the index of the electron
     *
     *  @return the index of the orbital that the electron occupies. For the bitset "100", this would be 2: the conversion from right-to-left is already made
     */
    size_t get_occupied_index(size_t electron_index) const { return occupation_indices(electron_index); }


    // PUBLIC METHODS
    /**
     *  Extracts the positions of the set bits from the this->unsigned_representation and places them in the this->occupation_indices
     */
    void updateOccupationIndices();

    /**
     *  @param p    the orbital index starting from 0, counted from right to left
     *
     *  @return if the p-th spatial orbital is occupied
     */
    bool isOccupied(size_t p) const;

    /**
     *  @param p    the orbital index starting from 0, counted from right to left
     *
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spatial orbital. Subsequently perform an in-place annihilation on the orbital p
     *
     *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
     */
    bool annihilate(size_t p);

    /**
     *  @param p        the orbital index starting from 0, counted from right to left
     *  @param sign     the current sign of the operator string
     *
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spatial orbital. Subsequently perform an in-place annihilation on the orbital p. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after annihilation.
     *
     *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
     */
    bool annihilate(size_t p, int& sign);

    /**
     *  @param indices      the indices of the orbitals that should be annihilated (the first index is annihilated first)
     *
     *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices p. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after the annihilations.
     *
     *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
     */
    bool annihilateAll(const std::vector<size_t>& indices, int& sign);

    /**
     *  @param p        the orbital index starting from 0, counted from right to left
     *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spatial orbital. Subsequently perform an in-place creation on the orbital p
     *
     *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
     */
    bool create(size_t p);

    /**
     *  @param p        the orbital index starting from 0, counted from right to left
     *  @param sign     the current sign of the operator string
     *
     *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spatial orbital. Subsequently perform an in-place creation on the orbital p. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after creation.
     *
     *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
     */
    bool create(size_t p, int& sign);

    /**
     *  @param p        the orbital index starting from 0, counted from right to left
     *
     *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on orbital p
     *
     *  Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
     */
    int operatorPhaseFactor(size_t p) const;

    /**
     *  @param index_start      the starting index (included), read from right to left
     *  @param index_end        the ending index (not included), read from right to left
     *
     *  @return the representation of a slice (i.e. a subset) of the spin string (read from right to left) between index_start (included) and index_end (not included)
     *
     *      Example:
     *          "010011".slice(1, 4) => "01[001]1" -> "001"
     */
    size_t slice(size_t index_start, size_t index_end) const;

    /**
     *  @param other        the other ONV
     *
     *  @return the number of different occupations between this ONV and the other, i.e. two times the number of electron excitations
     */
    size_t countNumberOfDifferences(const ONV& other) const;

    /**
     *  @param other        the other ONV
     *
     *  @return the indices of the orbitals (from right to left) that are occupied in this ONV, but unoccupied in the other
     */
    std::vector<size_t> findDifferentOccupations(const ONV &other) const;

    /**
     *  @param other        the other ONV
     *
     *  @return the indices of the orbitals (from right to left) that are occupied both this ONV and the other
     */
    std::vector<size_t> findMatchingOccupations(const ONV& other) const;

    /**
     * @return a string representation of the ONV
     */
    std::string asString() const;
};


}  // namespace GQCP

#endif  // GQCP_ONV_HPP
