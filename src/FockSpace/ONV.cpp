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
#include "FockSpace/ONV.hpp"

#include <boost/dynamic_bitset.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K                        the number of orbitals
 *  @param N                        the number of electrons
 *  @param unsigned_representation  the representation for the ONV as an unsigned integer
 */
ONV::ONV(size_t K, size_t N, size_t unsigned_representation) :
    ONV(K, N)
{
    this->unsigned_representation = unsigned_representation;
    this->updateOccupationIndices();  // throws error if the representation and N are not compatible
}


/**
 *  Constructs a default ONV without a representation
 *
 *  @param K                        the number of orbitals
 *  @param N                        the number of electrons
 */
ONV::ONV(size_t K, size_t N) :
    K (K),
    N (N),
    occupation_indices (VectorXs::Zero(N))
{
    this->occupation_indices = VectorXs::Zero(N);
}



/*
 *  OPERATORS
 */

/**
 *  @param os       the output stream which the ONV should be concatenated to
 *  @param onv      the ONV that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const ONV& onv) {
    return os<< onv.asString();
}


/**
 *  @param other    the other ONV
 *
 *  @return if this ONV is the same as the other ONV
 */
bool ONV::operator==(ONV& other) const {
    return this->unsigned_representation == other.unsigned_representation && this->K == other.K;  // this ensures that N, K and representation are equal
}


/**
 *  @param other    the other ONV
 *
 *  @return if this ONV is not the same as the other ONV
 */
bool ONV::operator!=(ONV& other) const {
    return !(this->operator==(other));
}



/*
 *  SETTERS
 */

/**
 *  @param unsigned_representation      the new representation as an unsigned integer
 *
 *  Set the representation of an ONV to a new representation and call update the occupation indices accordingly
 */
void ONV::set_representation(size_t unsigned_representation) {
    this->unsigned_representation = unsigned_representation;
    this->updateOccupationIndices();
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Extracts the positions of the set bits from the this->unsigned_representation and places them in the this->occupation_indices
 */
void ONV::updateOccupationIndices() {
    size_t l = this->unsigned_representation;
    int representation_electron = 0;
    while (l != 0) {
        this->occupation_indices(representation_electron) = __builtin_ctzl(l);  // retrieves occupation index
        representation_electron++;
        l ^= (l & -l);  // flip the least significant bit
    }
    if (representation_electron != this->N) {
        throw std::invalid_argument("The current representation and electron count are not compatible");
    }
}


/**
 *  @param p    the orbital index starting from 0, counted from right to left
 *
 *  @return if the p-th spatial orbital is occupied
 */
bool ONV::isOccupied(size_t p) const {

    if (p > this->K - 1) {
        throw std::invalid_argument("The index is out of the bitset bounds");
    }

    size_t operator_string = 1U << p;
    return this->unsigned_representation & operator_string;
}


/**
 *  @param indices      the orbital indices (starting from 0)
 *
 *  @return if all given indices are occupied
 */
bool ONV::areOccupied(const std::vector<size_t>& indices) const {

    // Check first if all indices are within bounds
    for (const auto& index : indices) {
        if (index > this->K - 1) {
            throw std::invalid_argument("The index is out of the bitset bounds");
        }
    }

    for (const auto& index : indices) {
        if (!this->isOccupied(index)) {
            return false;
        }
    }

    // Only if all indices have been tested to be occupied, we can return true
    return true;
}


/**
 *  @param p    the orbital index starting from 0, counted from right to left
 *
 *  @return if the p-th spatial orbital is not occupied
 */
bool ONV::isUnoccupied(size_t p) const {
    return !this->isOccupied(p);
}


/**
 *  @param indices      the orbital indices (starting from 0)
 *
 *  @return if all the given indices are unoccupied
 */
bool ONV::areUnoccupied(const std::vector<size_t>& indices) const {

    // Check first if all indices are within bounds
    for (const auto& index : indices) {
        if (index > this->K - 1) {
            throw std::invalid_argument("The index is out of the bitset bounds");
        }
    }

    for (const auto& index : indices) {
        if (this->isOccupied(index)) {
            return false;
        }
    }

    // Only if all indices have been tested to be unoccupied, we can return true
    return true;
}


/**
 *  @param index_start      the starting index (included), read from right to left
 *  @param index_end        the ending index (not included), read from right to left
 *
 *  @return the representation of a slice (i.e. a subset) of the spin string (read from right to left) between index_start (included) and index_end (not included)
 *
 *      Example:
 *          "010011".slice(1, 4) => "01[001]1" -> "001"
 */
size_t ONV::slice(size_t index_start, size_t index_end) const {

    // First, do some checks
    if (index_end <= index_start) {
        throw std::invalid_argument("index_end should be larger than index_start.");
    }

    if (index_end > this->K + 1) {
        throw std::invalid_argument("The last slicing index index_end cannot be greater than the number of spatial orbitals K.");
    }

    // The union of these conditions also include the case that index_start > this->K


    // Shift bits to the right
    size_t u = this->unsigned_representation >> index_start;


    // Create the correct mask
    size_t mask_length = index_end - index_start;
    size_t mask = ((1U) << mask_length) - 1;


    // Use the mask
    return u & mask;
}


/**
 *  @param p        the orbital index starting from 0, counted from right to left
 *
 *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on orbital p
 *
 *  Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
 */
int ONV::operatorPhaseFactor(size_t p) const {

    if (p == 0) {  // we can't give this to this->slice(0, 0)
        return 1;
    }
    size_t m = __builtin_popcountl(this->slice(0, p));  // count the number of set bits in the slice [0,p-1]

    if ( m % 2 == 0 ) {  // even number of electrons: phase factor (+1)
        return 1;
    } else {  // odd number of electrons: phase factor (-1)
        return -1;
    }
}


/**
 *  @param p    the orbital index starting from 0, counted from right to left
 *
 *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spatial orbital. Subsequently perform an in-place annihilation on the orbital p
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::annihilate(size_t p) {

    if (this->isOccupied(p)) {
        size_t operator_string = 1U << p;
        this->unsigned_representation &= ~operator_string;
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices      the orbital indices (starting from 0)
 *
 *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::annihilateAll(const std::vector<size_t>& indices) {

    if (this->areOccupied(indices)) {  // only if all indices are occupied, we will annihilate
        for (const auto& index : indices) {
            this->annihilate(index);
        }
        return true;
    } else {
        return false;
    }
}


/**
 *  @param p        the orbital index starting from 0, counted from right to left
 *  @param sign     the current sign of the operator string
 *
 *  @return if we can apply the annihilation operator (i.e. 1->0) for the p-th spatial orbital. Subsequently perform an in-place annihilation on the orbital p. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after annihilation.
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::annihilate(size_t p, int& sign) {

    if (this->annihilate(p)) {  // we have to first check if we can annihilate before applying the phase factor
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices      the indices of the orbitals that should be annihilated (the first index is annihilated first)
 *  @param sign     the current sign of the operator string
 *
 *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after the annihilations.
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::annihilateAll(const std::vector<size_t>& indices, int& sign) {

    if (this->areOccupied(indices)) {  // only if all indices are occupied, we will annihilate
        for (const auto& index : indices) {
            this->annihilate(index, sign);
        }
        return true;
    } else {
        return false;
    }
}


/**
 *  @param p        the orbital index starting from 0, counted from right to left
 *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spatial orbital. Subsequently perform an in-place creation on the orbital p
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::create(size_t p) {

    if (!this->isOccupied(p)) {
        size_t operator_string = 1U << p;
        this->unsigned_representation ^= operator_string;
        return true;
    } else {
        return false;
    }
}


/**
 *  @param p        the orbital index starting from 0, counted from right to left
 *  @param sign     the current sign of the operator string
 *
 *  @return if we can apply the creation operator (i.e. 0->1) for the p-th spatial orbital. Subsequently perform an in-place creation on the orbital p. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after creation.
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::create(size_t p, int& sign) {

    if (this->create(p)) {  // we have to first check if we can create before applying the phase factor
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices      the indices of the orbitals that should be created
 *
 *  @return if we can apply all creation operators (i.e. 0->1) on the given indices. Subsequently perform in-place creations on the given indices
 */
bool ONV::createAll(const std::vector<size_t>& indices) {

    if (this->areUnoccupied(indices)) {
        for (const auto& index : indices) {
            this->create(index);
        }
        return true;
    } else {
        return false;
    }
}

    
/**
 *  @param indices      the indices of the orbitals that should be annihilated (the first index is annihilated first)
 *  @param sign     the current sign of the operator string
 *
 *  @return if we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after the annihilations.
 *
 *  IMPORTANT: does not update the occupation indices for performance reasons, if required call updateOccupationIndices()!
 */
bool ONV::createAll(const std::vector<size_t>& indices, int& sign) {

    if (this->areUnoccupied(indices)) {
        for (const auto& index : indices) {
            this->create(index, sign);
        }
        return true;
    } else {
        return false;
    }
}


/**
 *  @param other        the other ONV
 *
 *  @return the number of different occupations between this ONV and the other, i.e. two times the number of electron excitations
 */
size_t ONV::countNumberOfDifferences(const ONV& other) const {
    return __builtin_popcountl(this->unsigned_representation ^ other.unsigned_representation);
}


/**
 *  @param other        the other ONV
 *
 *  @return the indices of the orbitals (from right to left) that are occupied in this ONV, but unoccupied in the other
 */
std::vector<size_t> ONV::findDifferentOccupations(const ONV &other) const {

    size_t differences = this->unsigned_representation ^ other.unsigned_representation;
    size_t occupied_differences = differences & this->unsigned_representation;  // this holds all indices occupied in this, but unoccupied in other

    size_t number_of_occupied_differences = __builtin_popcountl(occupied_differences);
    std::vector<size_t> positions (number_of_occupied_differences);


    // Find the positions of the set bits in occupied_differences
    for (size_t counter = 0; counter < number_of_occupied_differences; counter++) {  // counts the number of occupied differences we have already encountered
        size_t position = __builtin_ctzl(occupied_differences);  // count trailing zeros
        positions[counter] = position;

        occupied_differences ^= occupied_differences & -occupied_differences;  // annihilate the least significant set bit
    }

    return positions;
}


/**
 *  @param other        the other ONV
 *
 *  @return the indices of the orbitals (from right to left) that are occupied both this ONV and the other
 */
std::vector<size_t> ONV::findMatchingOccupations(const ONV& other) const {

    size_t matches = this->unsigned_representation & other.unsigned_representation;
    size_t number_of_occupied_matches = __builtin_popcountl(matches);
    Vectoru positions (number_of_occupied_matches);


    // Find the positions of the set bits in occupied_differences
    for (size_t counter = 0; counter < number_of_occupied_matches; counter++) {  // counts the number of occupied differences we have already encountered
        size_t position = __builtin_ctzl(matches);  // count trailing zeros
        positions[counter] = position;

        matches ^= matches & -matches;  // annihilate the least significant set bit
    }

    return positions;
}


/**
 *  @return a string representation of the ONV
 */
std::string ONV::asString() const {
    boost::dynamic_bitset<> transfer_set (this->K, this->unsigned_representation);
    std::string buffer;
    boost::to_string(transfer_set, buffer);
    return buffer;
}


}  // namespace GQCP
