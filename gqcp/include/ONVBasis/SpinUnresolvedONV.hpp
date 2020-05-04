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


#include "Mathematical/Representation/Matrix.hpp"


namespace GQCP {


/**
 *  A spin-unresolved occupation number vector.

 *  An spin-unresolved ONV in quantum chemistry is a string of spinor creation operators acting on top of a vacuum state.
 *  An example for 3 electrons in a general spinor basis spanned by 4 spinors is
 *      a_1^\dagger a_2^\dagger a_3^\dagger |vac> = |1,1,1,0>
 *
 *  In GQCP, bitstrings are read from right to left. This means that the least significant bit relates to the first orbital.
 *  This notation is consistent with how normally bit strings are read in computer science, leading to more efficient code. The least significant bit has index 0. The previous example is then represented by the bit string "0111" (with unsigned representation 7).
 * 
 *  IMPORTANT: Since this type heavily uses unsigned bitset representations, it should not be used for M > 64 spinors.
 */
class SpinUnresolvedONV {
private:
    size_t M;  // the number of spinors that this ONV is expressed in
    size_t N;  // the number of electrons that appear in this ONV, i.e. the number of spinor that is occupied

    size_t unsigned_representation;  // the representation of this ONV as an unsigned integer
    VectorXs occupation_indices;     // the indices of the spinors that are occupied in this ONV
                                     // it is a vector of N elements in which occupation_indices[j] returns the index of the spinor that the electron j occupies


public:
    // CONSTRUCTORS

    /**
     *  Create a SpinResolvedONV from an unsigned representation.
     * 
     *  @param M                                the number of spinors that this ONV is expressed in
     *  @param N                                the number of electrons that appear in this ONV, i.e. the number of spinor that is occupied
     *  @param unsigned_representation          the representation of this ONV as an unsigned integer
     */
    SpinUnresolvedONV(const size_t M, const size_t N, const size_t unsigned_representation);

    /**
     *  Construct a SpinResolvedONV ONV without an unsigned representation
     *
     *  @param M                the number of spinors that this ONV is expressed in
     *  @param N                the number of electrons that appear in this ONV, i.e. the number of spinor that is occupied
     */
    SpinUnresolvedONV(const size_t M, const size_t N);


    // OPERATORS

    /**
     *  @param os       the output stream which the spin-unresolved ONV should be concatenated to
     *  @param onv      the spin-unresolved ONV that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const SpinUnresolvedONV& onv);

    /**
     *  @param other    the other spin-unresolved ONV
     *
     *  @return if this spin-unresolved ONV is the same as the other spin-unresolved ONV
     */
    bool operator==(const SpinUnresolvedONV& other) const;

    /**
     *  @param other    the other spin-unresolved ONV
     *
     *  @return if this spin-unresolved ONV is not the same as the other spin-unresolved ONV
     */
    bool operator!=(const SpinUnresolvedONV& other) const;


    // SETTERS
    /**
     *  @param unsigned_representation      the new representation as an unsigned integer
     *
     *  Set the representation of an spin-unresolved ONV to a new representation and call update the occupation indices accordingly
     */
    void set_representation(const size_t unsigned_representation);


    // GETTERS
    size_t get_unsigned_representation() const { return unsigned_representation; }
    const VectorXs& get_occupation_indices() const { return occupation_indices; }

    /**
     *  @param electron_index       the index of the electron
     *
     *  @return the index of the orbital that the electron occupies. For the bitset "100", this would be 2: the conversion from right-to-left is already made
     */
    size_t get_occupation_index(const size_t electron_index) const { return occupationIndexOf(electron_index); }


    // PUBLIC METHODS

    /**
     *  Annihilate the electron at a given spinor index.
     * 
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spinor and subsequently perform an in-place annihilation on that spinor
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool annihilate(const size_t p);

    /**
     *  Annihilate the electron at a given spinor index, keeping track of any sign changes.
     * 
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *  @param sign         the current sign of the operator string
     *
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spinor and subsequently perform an in-place annihilation on that spinor. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after annihilation
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool annihilate(const size_t p, int& sign);

    /**
     *  Annihilate the electrons at the given spinor indices.
     * 
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *
     *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. If possible, subsequently perform in-place annihilations on all the given indices.
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool annihilateAll(const std::vector<size_t>& indices);

    /**
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *  @param sign             the current sign of the operator string
     *
     *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the operator string after the annihilations.
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool annihilateAll(const std::vector<size_t>& indices, int& sign);

    /**
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *
     *  @return if all the spinors with the given indices are occupied
     */
    bool areOccupied(const std::vector<size_t>& indices) const;

    /**
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *
     *  @return if all the spinors with the given indices are unoccupied
     */
    bool areUnoccupied(const std::vector<size_t>& indices) const;

    /**
     *  @return a string representation of this spin-unresolved ONV
     */
    std::string asString() const;

    /**
     *  @param other        the other spin-unresolved ONV
     *
     *  @return the number of different occupations between this spin-unresolved ONV and the other, i.e. two times the number of electron excitations
     */
    size_t countNumberOfDifferences(const SpinUnresolvedONV& other) const;

    /**
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     * 
     *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spinor and subsequently perform an in-place annihilation on that spinor
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool create(size_t p);

    /**
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *  @param sign         the current sign of the operator string
     *
     *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spinor and subsequently perform an in-place annihilation on that spinor. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after creation.
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool create(size_t p, int& sign);

    /**
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *
     *  @return if we can apply all creation operators (i.e. 0->1) on the given indices. Subsequently perform in-place creations on the given indices
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool createAll(const std::vector<size_t>& indices);

    /**
     *  @param indices          the 0-based spinor indices, counted in this ONV from right to left
     *  @param sign             the current sign of the operator string
     *
     *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the operator string after the annihilations.
     *
     *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
     */
    bool createAll(const std::vector<size_t>& indices, int& sign);

    /**
     *  @param other            the other spin-unresolved ONV
     *
     *  @return the indices of the spinors (from right to left) that are occupied in this spin-unresolved ONV, but unoccupied in the other
     */
    std::vector<size_t> findDifferentOccupations(const SpinUnresolvedONV& other) const;

    /**
     *  @param other            the other spin-unresolved ONV
     *
     *  @return the indices of the spinors (from right to left) that are occupied both this spin-unresolved ONV and the other
     */
    std::vector<size_t> findMatchingOccupations(const SpinUnresolvedONV& other) const;

    /**
     *  Iterate over every occupied spinor index in this ONV and apply the given callback function.
     * 
     *  @param callback         the function that should be called in every iteration step over all occupied spinor indices. The argument of this callback function is the index of the occupied spinor.
     */
    void forEach(const std::function<void(const size_t)>& callback) const;

    /**
     *  Iterate over every unique pair of occupied spinor indices in this ONV and apply the given callback function.
     * 
     *  @param callback         the function that should be called in every iteration step over all pairs of occupied spinor indices. The arguments of this callback function are the indices of the occupied spinor, where the first index is always larger than the second.
     */
    void forEach(const std::function<void(const size_t, const size_t)>& callback) const;

    /**
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *
     *  @return if the p-th spinor is occupied
     */
    bool isOccupied(const size_t p) const;

    /**
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *
     *  @return if the p-th spinor is not occupied
     */
    bool isUnoccupied(const size_t p) const;

    /**
     *  @return the number of electrons that this ONV contains.
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @param electron_index           the index of an electron in this ONV
     *  
     *  @return the index of spinor that the given electron (index) occupies
     */
    size_t occupationIndexOf(const size_t electron_index) const { return occupation_indices(electron_index); }

    /**
     *  @param p            the 0-based spinor index, counted in this ONV from right to left
     *
     *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on spinor p
     *
     *  @example Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
     */
    int operatorPhaseFactor(const size_t p) const;

    /**
     *  @param index_start      the starting index (included), read from right to left
     *  @param index_end        the ending index (not included), read from right to left
     *
     *  @return the representation of a slice (i.e. a subset) of the ONV (read from right to left) between index_start (included) and index_end (not included)
     *
     *  @example
     *      "010011".slice(1, 4) => "01[001]1" -> "001", where the spinor indices are:
     *       543210
     */
    size_t slice(const size_t index_start, const size_t index_end) const;

    /**
     *  Extract the positions of the set bits from 'this->unsigned_representation' and places them in 'this->occupation_indices'.
     */
    void updateOccupationIndices();
};


}  // namespace GQCP
