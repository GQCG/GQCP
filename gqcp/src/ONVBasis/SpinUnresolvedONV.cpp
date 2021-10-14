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

#include "ONVBasis/SpinUnresolvedONV.hpp"

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"

#include <boost/dynamic_bitset.hpp>

#include <algorithm>
#include <numeric>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param M                                The number of spinors that this ONV is expressed in.
 *  @param N                                The number of electrons that appear in this ONV, i.e. the number of spinors that is occupied.
 *  @param unsigned_representation          The representation of this ONV as an unsigned integer.
 */
SpinUnresolvedONV::SpinUnresolvedONV(const size_t M, const size_t N, const size_t unsigned_representation) :
    SpinUnresolvedONV(M, N) {

    this->unsigned_representation = unsigned_representation;
    this->updateOccupationIndices();  // Throws error if the representation and N are not compatible.
}


/**
 *  Construct a SpinResolvedONV ONV without an unsigned representation.
 *
 *  @param M                The number of spinors that this ONV is expressed in.
 *  @param N                The number of electrons that appear in this ONV, i.e. the number of spinors that is occupied.
 */
SpinUnresolvedONV::SpinUnresolvedONV(const size_t M, const size_t N) :
    M {M},
    N {N},
    occupied_indices(N, 0) {}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Create a spin-unresolved ONV from a textual/string representation.
 *
 *  @param string_representation                The textual representation of the spin-unresolved ONV, for example "0011", indicating that the first two spinors should be occupied.
 *
 *  @return A spin-unresolved ONV from a textual/string representation.
 */
SpinUnresolvedONV SpinUnresolvedONV::FromString(const std::string& string_representation) {

    boost::dynamic_bitset<> intermediate_bitset {string_representation};  // The least significant bit has the highest position in the string.
    const auto M = intermediate_bitset.size();
    const auto N = intermediate_bitset.count();

    const auto unsigned_representation = static_cast<size_t>(intermediate_bitset.to_ulong());

    return SpinUnresolvedONV(M, N, unsigned_representation);
}


/**
 *  Create a spin-unresolved ONV from a set of occupied indices.
 *
 *  @param occupied_indices             The indices that the electrons occupy, in order: e.g. the i-th element describes the spinor that the i-th electron occupies.
 *  @param M                            The total number of spinors.
 *
 *  @return A spin-resolved ONV from a set of occupied indices.
 */
SpinUnresolvedONV SpinUnresolvedONV::FromOccupiedIndices(const std::vector<size_t>& occupied_indices, const size_t M) {

    // Generate the corresponding unsigned representation and use that constructor.
    size_t unsigned_representation = 0;
    for (const auto& index : occupied_indices) {
        unsigned_representation += std::pow(2, index);
    }

    const size_t N = occupied_indices.size();
    return SpinUnresolvedONV(M, N, unsigned_representation);
}


/**
 *  Create a spin-unresolved ONV that represents the GHF single Slater determinant, occupying the N spinors with the lowest spinor energy.
 *
 *  @param M                            The number of spinors.
 *  @param N                            The number of electrons.
 *  @param orbital_energies             The single-particle energies of the spinors.
 *
 *  @return A spin-unresolved ONV that represents the GHF single Slater determinant.
 */
SpinUnresolvedONV SpinUnresolvedONV::GHF(const size_t M, const size_t N, const VectorX<double>& orbital_energies) {

    // The GHF ONV is that one in which the N spinors with the lowest energy are occupied.

    // Create an array that contains the indices of the spinors with ascending energy.
    std::vector<size_t> indices(M);                // Zero-initialized with M elements.
    std::iota(indices.begin(), indices.end(), 0);  // Start with 0.

    // Sort the indices according to the orbital energies.
    std::stable_sort(indices.begin(), indices.end(), [&orbital_energies](const size_t i, const size_t j) { return orbital_energies(i) < orbital_energies(j); });

    const std::vector<size_t> occupied_indices {indices.begin(), indices.begin() + N};  // The first N elements.
    return SpinUnresolvedONV::FromOccupiedIndices(occupied_indices, M);
}


/*
 *  OPERATORS
 */

/**
 *  @param os       The output stream which the ONV should be concatenated to.
 *  @param onv      The ONV that should be concatenated to the output stream.
 *
 *  @return The updated output stream.
 */
std::ostream& operator<<(std::ostream& os, const SpinUnresolvedONV& onv) {
    return os << onv.asString();
}


/**
 *  @param other    The other ONV.
 *
 *  @return If this ONV is the same as the other ONV.
 */
bool SpinUnresolvedONV::operator==(const SpinUnresolvedONV& other) const {
    return this->unsigned_representation == other.unsigned_representation && this->M == other.M;  // this ensures that N, M and representation are equal
}


/**
 *  @param other    The other ONV.
 *
 *  @return If this ONV is not the same as the other ONV.
 */
bool SpinUnresolvedONV::operator!=(const SpinUnresolvedONV& other) const {
    return !(this->operator==(other));
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Annihilate the electron at a given spinor index.
 *
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *
 *  @return If we can apply the annihilation operator (i.e. 1->0) for the p-th spinor and subsequently perform an in-place annihilation on that spinor.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::annihilate(const size_t p) {

    if (this->isOccupied(p)) {
        const size_t operator_string = 1U << p;
        this->unsigned_representation &= ~operator_string;
        return true;
    } else {
        return false;
    }
}


/**
 *  Annihilate the electron at a given spinor index, keeping track of any sign changes.
 *
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *  @param sign         The current sign of the operator string.
 *
 *  @return If we can apply the annihilation operator (i.e. 1->0) for the p-th spinor and subsequently perform an in-place annihilation on that spinor. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after annihilation.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::annihilate(const size_t p, int& sign) {

    if (this->annihilate(p)) {  // We have to first check if we can annihilate before applying the phase factor.
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}


/**
 *  Annihilate the electrons at the given spinor indices.
 *
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *
 *  @return If we can apply all annihilation operators (i.e. 1->0) on the given indices. If possible, subsequently perform in-place annihilations on all the given indices.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::annihilateAll(const std::vector<size_t>& indices) {

    if (this->areOccupied(indices)) {  // Only if all indices are occupied, we will annihilate.
        for (const auto& index : indices) {
            this->annihilate(index);
        }
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *  @param sign             The current sign of the operator string.
 *
 *  @return If we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the operator string after the annihilations.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::annihilateAll(const std::vector<size_t>& indices, int& sign) {

    if (this->areOccupied(indices)) {  // Only if all indices are occupied, we will annihilate.
        for (const auto& index : indices) {
            this->annihilate(index, sign);
        }
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *
 *  @return If all the spinors with the given indices are occupied.
 */
bool SpinUnresolvedONV::areOccupied(const std::vector<size_t>& indices) const {

    for (const auto& index : indices) {
        if (!this->isOccupied(index)) {
            return false;
        }
    }

    // Only if all indices have been tested to be occupied, we can return true.
    return true;
}


/**
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *
 *  @return If all the spinors with the given indices are unoccupied.
 */
bool SpinUnresolvedONV::areUnoccupied(const std::vector<size_t>& indices) const {

    for (const auto& index : indices) {
        if (this->isOccupied(index)) {
            return false;
        }
    }

    // Only if all indices have been tested to be unoccupied, we can return true.
    return true;
}


/**
 *  @return A string representation of this spin-unresolved ONV.
 */
std::string SpinUnresolvedONV::asString() const {

    boost::dynamic_bitset<> intermediate_bitset {this->M, this->unsigned_representation};
    std::string text;  // The string that will contain the textual representation of this spin-unresolved ONV.

    boost::to_string(intermediate_bitset, text);
    std::reverse(text.begin(), text.end());  // The bitstring is read from right to left, but we return it from left to right for easier reading.

    return text;
}


/**
 *  Calculate the overlap <on|of>: the projection of between this spin-unresolved ONV ('of') and another spin-unresolved ONV ('on'), expressed in different general orthonormal spinor bases.
 *
 *  @param onv_on                       The spin-unresolved ONV that should be projected on.
 *  @param C_of                         The transformation between the spinors related to the ONV that is being projected and the atomic spinors.
 *  @param C_on                         The transformation between the spinors related to the ONV that is being projected on and the atomic spinors.
 *  @param S                            The overlap matrix of the underlying AOs.
 *
 *  @return The overlap element <on|of>.
 */
double SpinUnresolvedONV::calculateProjection(const SpinUnresolvedONV& onv_on, const GTransformation<double>& C_of, const GTransformation<double>& C_on, const SquareMatrix<double>& S) const {

    // Make a reference copy in order to improve readibility of the following code.
    const auto& onv_of = *this;


    // Calculate the raw matrix representation of the transformation between the spinor bases, since we're going to have to slice its rows and columns.
    MatrixX<double> U = C_on.matrix().adjoint() * S * C_of.matrix();


    // U's columns should be the ones occupied in the 'of'-ONV.
    // U's rows should be the ones occupied in the 'on'-ONV.
    // While waiting for Eigen 3.4 to release (which has better slicing APIs), we'll remove the UNoccupied rows/columns.
    const auto unoccupied_indices_of = onv_of.unoccupiedIndices();
    const auto unoccupied_indices_on = onv_on.unoccupiedIndices();

    U.removeColumns(unoccupied_indices_of);
    U.removeRows(unoccupied_indices_on);


    // The requested overlap element is the determinant of the resulting matrix.
    return U.determinant();
}


/**
 *  @param other        The other ONV.
 *
 *  @return The number of different occupations between this ONV and the other.
 */
size_t SpinUnresolvedONV::countNumberOfDifferences(const SpinUnresolvedONV& other) const {
    return __builtin_popcountl(this->unsigned_representation ^ other.unsigned_representation);
}


/**
 *  @param other        The other ONV.
 *
 *  @return The number of electron excitations between this ONV and the other.
 */
size_t SpinUnresolvedONV::countNumberOfExcitations(const SpinUnresolvedONV& other) const {

    const auto differences = this->unsigned_representation ^ other.unsigned_representation;
    return std::min(__builtin_popcountl(this->unsigned_representation & differences), __builtin_popcountl(other.unsigned_representation & differences));
}


/**
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *
 *  @return If we can apply the creation operator (i.e. 0->1) for the p-th spinor and subsequently perform an in-place annihilation on that spinor.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::create(const size_t p) {

    if (!this->isOccupied(p)) {
        const size_t operator_string = 1U << p;
        this->unsigned_representation ^= operator_string;
        return true;
    } else {
        return false;
    }
}


/**
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *  @param sign         The current sign of the operator string.
 *
 *  @return If we can apply the creation operator (i.e. 0->1) for the p-th spinor and subsequently perform an in-place annihilation on that spinor. Furthermore, update the sign according to the sign change (+1 or -1) of the spin string after creation.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::create(const size_t p, int& sign) {

    if (this->create(p)) {  // We have to first check if we can create before applying the phase factor.
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}


/**
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *
 *  @return If we can apply all creation operators (i.e. 0->1) on the given indices. Subsequently perform in-place creations on the given indices.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::createAll(const std::vector<size_t>& indices) {

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
 *  @param indices          The 0-based spinor indices, counted in this ONV from right to left.
 *  @param sign             The current sign of the operator string.
 *
 *  @return If we can apply all annihilation operators (i.e. 1->0) on the given indices. Subsequently perform in-place annihilations on the given indices. Furthermore, update the sign according to the sign change (+1 or -1) of the operator string after the annihilations.
 *
 *  IMPORTANT: This function does not update the occupation indices for performance reasons. If required, call updateOccupationIndices()!
 */
bool SpinUnresolvedONV::createAll(const std::vector<size_t>& indices, int& sign) {

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
 *  @param other            The other spin-unresolved ONV.
 *
 *  @return The indices of the spinors (from right to left) that are occupied in this spin-unresolved ONV, but unoccupied in the other.
 */
std::vector<size_t> SpinUnresolvedONV::findDifferentOccupations(const SpinUnresolvedONV& other) const {

    const size_t differences = this->unsigned_representation ^ other.unsigned_representation;
    size_t occupied_differences = differences & this->unsigned_representation;  // This holds all indices occupied in this, but unoccupied in other.

    const size_t number_of_occupied_differences = __builtin_popcountl(occupied_differences);
    std::vector<size_t> positions;
    positions.reserve(number_of_occupied_differences);

    // Find the positions of the set bits in occupied_differences.
    for (size_t counter = 0; counter < number_of_occupied_differences; counter++) {  // Counts the number of occupied differences we have already encountered.
        const size_t position = __builtin_ctzl(occupied_differences);                // Count trailing zeros.
        positions.push_back(position);

        occupied_differences ^= occupied_differences & -occupied_differences;  // Annihilate the least significant set bit.
    }

    return positions;
}


/**
 *  @param other            The other spin-unresolved ONV.
 *
 *  @return The indices of the spinors (from right to left) that are occupied both this spin-unresolved ONV and the other.
 */
std::vector<size_t> SpinUnresolvedONV::findMatchingOccupations(const SpinUnresolvedONV& other) const {

    size_t matches = this->unsigned_representation & other.unsigned_representation;

    const size_t number_of_occupied_matches = __builtin_popcountl(matches);
    std::vector<size_t> positions;
    positions.reserve(number_of_occupied_matches);


    // Find the positions of the set bits in occupied_differences.
    for (size_t counter = 0; counter < number_of_occupied_matches; counter++) {  // Counts the number of occupied differences we have already encountered.
        const size_t position = __builtin_ctzl(matches);                         // Count trailing zeros.
        positions.push_back(position);

        matches ^= matches & -matches;  // Annihilate the least significant set bit.
    }

    return positions;
}


/**
 *  Iterate over every occupied spinor index in this ONV and apply the given callback function.
 *
 *  @param callback         The function that should be called in every iteration step over all occupied spinor indices. The argument of this callback function is the index of the occupied spinor.
 */
void SpinUnresolvedONV::forEach(const std::function<void(const size_t)>& callback) const {

    // Loop over every electron in this ONV and retrieve the index of the spinor that it occupies.
    for (size_t e = 0; e < this->numberOfElectrons(); e++) {
        const size_t p = this->occupationIndexOf(e);
        callback(p);
    }  // electron index loop
}


/**
 *  Iterate over every unique pair of occupied spinor indices in this ONV and apply the given callback function.
 *
 *  @param callback         The function that should be called in every iteration step over all pairs of occupied spinor indices. The arguments of this callback function are the indices of the occupied spinor, where the first index is always larger than the second.
 */
void SpinUnresolvedONV::forEach(const std::function<void(const size_t, const size_t)>& callback) const {

    // Loop over every electron in this ONV and retrieve the index of the spinor that it occupies.
    for (size_t e1 = 0; e1 < this->numberOfElectrons(); e1++) {
        const size_t p = this->occupationIndexOf(e1);

        // Loop over every different electron in this ONV and retrieve the index of the spinor that it occupies.
        for (size_t e2 = 0; e2 < e1; e2++) {
            const size_t q = this->occupationIndexOf(e2);

            callback(p, q);
        }  // electron 2 index loop
    }      // electron 1 index loop
}


/**
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *
 *  @return If the p-th spinor is occupied.
 */
bool SpinUnresolvedONV::isOccupied(const size_t p) const {

    if (p > this->M - 1) {
        throw std::invalid_argument("SpinUnresolvedONV::isOccupied(size_t): The index is out of the bitset bounds.");
    }

    const size_t operator_string = 1U << p;
    return this->unsigned_representation & operator_string;
}


/**
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *
 *  @return If the p-th spinor is not occupied.
 */
bool SpinUnresolvedONV::isUnoccupied(const size_t p) const {
    return !this->isOccupied(p);
}


/**
 *  @param p            The 0-based spinor index, counted in this ONV from right to left.
 *
 *  @return The phase factor (+1 or -1) that arises by applying an annihilation or creation operator on spinor p.
 *
 *  @example Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
 */
int SpinUnresolvedONV::operatorPhaseFactor(const size_t p) const {

    if (p == 0) {  // We can't give this to this->slice(0, 0).
        return 1;
    }
    const size_t m = __builtin_popcountl(this->slice(0, p));  // Count the number of set bits in the slice [0,p-1].

    if (m % 2 == 0) {  // Even number of electrons: phase factor (+1).
        return 1;
    } else {  // Odd number of electrons: phase factor (-1).
        return -1;
    }
}


/**
 *  @return The implicit orbital space that is related to this spin-unresolved ONV by taking this as a reference determinant.
 */
OrbitalSpace SpinUnresolvedONV::orbitalSpace() const {

    // Create an occupied-virtual orbital space.
    return OrbitalSpace(this->occupiedIndices(), this->unoccupiedIndices());
}


/**
 *  @param unsigned_representation      The new representation as an unsigned integer.
 *
 *  Set the representation of an ONV to a new representation and call update the occupation indices accordingly.
 */
void SpinUnresolvedONV::replaceRepresentationWith(const size_t unsigned_representation) {

    this->unsigned_representation = unsigned_representation;
    this->updateOccupationIndices();
}


/**
 *  @param index_start      The starting index (included), read from right to left.
 *  @param index_end        The ending index (not included), read from right to left.
 *
 *  @return The representation of a slice (i.e. a subset) of the ONV (read from right to left) between index_start (included) and index_end (not included).
 *
 *  @example
 *      "010011".slice(1, 4) => "01[001]1" -> "001", where the spinor indices are:
 *       543210
 */
size_t SpinUnresolvedONV::slice(const size_t index_start, const size_t index_end) const {

    // First, do some checks
    if (index_end <= index_start) {
        throw std::invalid_argument("SpinUnresolvedONV::slice(size_t, size_t): index_end should be larger than index_start.");
    }

    if (index_end > this->M + 1) {
        throw std::invalid_argument("SpinUnresolvedONV::slice(size_t, size_t): The last slicing index index_end cannot be greater than the number of spatial orbitals M.");
    }

    // The union of these conditions also include the case that index_start > this->M.


    // Shift bits to the right.
    const size_t u = this->unsigned_representation >> index_start;


    // Create the correct mask.
    const size_t mask_length = index_end - index_start;
    const size_t mask = ((1U) << mask_length) - 1;


    // Use the mask.
    return u & mask;
}


/**
 *  @return The spinor indices that are not occupied in this ONV.
 */
std::vector<size_t> SpinUnresolvedONV::unoccupiedIndices() const {

    // Create a vector containing all indices.
    std::vector<size_t> all_indices(this->M);
    std::iota(all_indices.begin(), all_indices.end(), 0);  // Fill all_indices with increasing numbers, starting by 0.

    // The unoccupied indices are {all indices}\{occupied indices}.
    std::vector<size_t> unoccupied_indices;
    std::set_difference(all_indices.begin(), all_indices.end(), this->occupied_indices.begin(), this->occupied_indices.end(),
                        std::inserter(unoccupied_indices, unoccupied_indices.begin()));

    return unoccupied_indices;
}


/**
 *  Extract the positions of the set bits from 'this->unsigned_representation' and places them in 'this->occupied_indices'.
 */
void SpinUnresolvedONV::updateOccupationIndices() {

    size_t representation_copy = this->unsigned_representation;
    int electron_index = 0;

    while (representation_copy != 0) {
        this->occupied_indices[electron_index] = __builtin_ctzl(representation_copy);  // Retrieves occupation index.
        electron_index++;

        representation_copy ^= (representation_copy & -representation_copy);  // Flip the least significant bit.
    }

    if (electron_index != this->N) {
        throw std::invalid_argument("SpinUnresolvedONV::updateOccupationIndices(): The current representation and electron count are not compatible.");
    }
}


}  // namespace GQCP
