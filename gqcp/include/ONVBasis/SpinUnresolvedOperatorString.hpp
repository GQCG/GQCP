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

#include <iterator>
#include <vector>


namespace GQCP {


/**
 *  A spin-unresolved operator string.

 *  A spin unresolved operator string represents a string of either annihilation or creation operators by its indices and a corresponding phase factor p.
 *  For example, an operator string represented by indices <1, 2, 3> represents either:
 *       p * a_1^\dagger a_2^\dagger a_3^\dagger
 *  or
 *       p * a_1 a_2 a_3.
 *  
 *  An operator string is always represented by the indices of the operators. Whether it denotes annihilation or creation operators depends on the context in which an operator string is used.
 *  The operator strings are different from ONV's. Instead of representing the way orbitals are occupied, they purely represent the order of certain operators.
 */
class SpinUnresolvedOperatorString {
private:
    // The vector representing the indices of the annihilation or creation operators in the operator string.
    std::vector<size_t> indices;

    // The phase factor associated with this string of annihilation or creation operators.
    int p;

public:
    /*
     * MARK: Constructors
     */

    /**
     *  Construct a `SpinUnresolvedOperatorString` from the vector of indices that it encapsulates.
     * 
     *  @param index_vector                The vector containing the operator indices.
     */
    SpinUnresolvedOperatorString(const std::vector<size_t>& index_vector) :
        indices {index_vector},
        p {1} {};


    /**
     *  Construct a `SpinUnresolvedOperatorString` from the vector of indices that it encapsulates and a specified phase factor.
     * 
     *  @param index_vector                The vector containing the operator indices.
     *  @param phase_factor                The phase factor associated with the operator string, in case fermion anti-commutation rules have been used to modify the operator string. 
     */
    SpinUnresolvedOperatorString(const std::vector<size_t>& index_vector, const int phase_factor) :
        indices {index_vector},
        p {phase_factor} {};


    /*
     * MARK: Named constructors
     */

    /**
     *  Construct a `SpinUnresolvedOperatorString` from the occupied indices of a `SpinUnresolvedONV`.
     * 
     *  @param onv                The `SpinUnresolvedONV` encapsulating the operator indices of the ONV.
     */
    static SpinUnresolvedOperatorString FromONV(const SpinUnresolvedONV& onv) { return SpinUnresolvedOperatorString(onv.occupiedIndices()); };


    /*
     *  MARK: General information
     */

    /**
     *  Retrieve the operator indices from the `SpinUnresolvedOperatorString`.
     * 
     *  @return The vector containing the operator indices.
     */
    const std::vector<size_t>& operatorIndices() const { return this->indices; }

    /**
     *  Retrieve the phase factor corresponding to the `SpinUnresolvedOperatorString`.
     * 
     *  @return The phase factor associated with the operator string.
     */
    int phaseFactor() const { return this->p; }

    /**
     *  Check whether the operator string in question will result in zero when applied to the wave function.
     *
     *  Note: There are different cases when an operator string will result in a zero value. This method checks all of them.
     * 
     *  @return Whether the operator string in question will result in zero when applied to the wave function.
     */
    bool isZero() const;


    /*
     *  MARK: Public methods
     */

    /**
     *  Retrieve the phase factor after sorting the `SpinUnresolvedOperatorString`.
     * 
     *  Note: Please refer to method `sort()` to perform the actual sorting.
     * 
     *  @return The phase factor after sorting the operator string.
     */
    int phaseFactorAfterSorting();

    /**
     *  Sort the operator string (in-place) in ascending order and adjust its phase factor. If two indices are equal, the phase factor may not be correct but it does not really matter since its effect on a wavefunction is zero.
     */
    void sort();

    /**
     *  Partition the `SpinUnresolvedOperatorString` into two new operator strings: a system and an environment.
     * 
     *  @param partition    The partition of the operator string into a system (denoted by 'I') and an environment (denoted by 'J').
     *  
     *  For example: Operator string "a1a2a4a0a3" is partitioned into system "a1a4a0" and environment "a2a3" by the partition {'I', 'J', 'I', 'I', 'J'}.
     */
    std::vector<SpinUnresolvedOperatorString> partitionIntoTwoSubsystems(const std::vector<char>& partition);
};


}  // namespace GQCP
