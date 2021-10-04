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


#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedOperatorString.hpp"
#include "QuantumChemical/Spin.hpp"
#include "QuantumChemical/SpinResolved.hpp"

#include <vector>


namespace GQCP {


/**
 *  A spin resolved operator string.
  *
 *  A spin resolved operator string represents a string of either annihilation or creation operators by their orbital indices. The indices can belong to either alpha or beta operators.
 *  For example, an operator string represented by index-spin pairs <1a, 1b, 2a, 3b> represents either:
 *      a_1_alpha^\dagger a_1_beta^\dagger a_2_alpha^\dagger a_3_beta^\dagger
 *  or
 *      a_1_alpha a_1_beta a_2_alpha a_3_beta.
 *  
 *  An operator string is always represented by pairs, containing the indices of the operators on the one hand and the spin (alpha or beta) of the operator on the other. Whether it denotes annihilation or creation operators depends on the context in which an operator string is used.
 *  The operator strings are different from ONV's. Instead of representing the way orbitals are occupied, they purely represent the order of certain operators.
 */
class SpinResolvedOperatorString {
private:
    // The vector representing the indices of the annihilation or creation operators in the operator string. Each index is paired with a certain spin.
    std::vector<std::pair<size_t, Spin>> index_spin_pairs;

public:
    /*
     * MARK: Constructors
     */

    /**
     *  Construct a `SpinResolvedOperatorString` from the vector of indices that it encapsulates, together with a vector containing the spins of the individual operators.
     * 
     *  @param index_vector                The vector containing the operator indices.
     *  @param spin_vector                 The vector containing the spins associated with each individual operator index.
     * 
     * Note: The index vector element zero will be matched with the spin vector element zero, the index vector element one will be matched with spin vector element one and so on.
     */
    SpinResolvedOperatorString(const std::vector<size_t>& index_vector, const std::vector<Spin>& spin_vector);


    /*
     * MARK: Named constructors
     */

    static SpinResolvedOperatorString FromONV(const SpinResolvedONV& onv);


    /*
     *  MARK: General information
     */

    /**
     *  Retrieve the operator indices from the `SpinResolvedOperatorString`.
     */
    std::vector<size_t> operatorIndices() const;


    /**
     *  Retrieve the operator spins from the `SpinResolvedOperatorString`.
     */
    std::vector<Spin> operatorSpins() const;

    /**
     *  Retrieve the number of operators in the `SpinResolvedOperatorString`.
     */
    size_t size() const { return this->index_spin_pairs.size(); }


    /*
     *  MARK: Spin component acces
     */

    /**
     *  Return the spin resolved operator string split in its two separate components, one containing the alpha operators, one containing the beta operators. The alpha operator string will recieve the total phase factor assiciated with this action, the phase factor of the beta string will be one.
     */
    SpinResolved<SpinUnresolvedOperatorString> spinResolve() const;
};


}  // namespace GQCP
