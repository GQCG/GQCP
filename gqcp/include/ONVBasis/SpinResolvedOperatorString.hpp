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


#include "ONVBasis/SpinUnresolvedOperatorString.hpp"
#include "QuantumChemical/Spin.hpp"
#include "QuantumChemical/SpinResolved.hpp"


namespace GQCP {


/**
 *  A spin resolved operator string.

 *  An spin resolved operator string represents a string of either annihilation or creation operators by its indices. The indices can belong to either alpha or beta operators.
 *  For example, an operator string represented by indices <1a, 1b, 2a, 3b> represents either:
 *      a_1_alpha^\dagger a_1_beta^\dagger a_2_alpha^\dagger a_3_beta^\dagger
 *  or
 *      a_1_alpha a_1_beta a_2_alpha a_3_beta.
 *  
 *  An operator string is always represented by pairs, containing the indices of the operators on the one hand and the spin (alpha or beta) on the other. Whether it denotes annihilation or creation operators depends on the context in which an operator string is used.
 *  The operator strings are different from ONV's. Instead of representing the way orbitals are occupied, they purely represent the order of certain operators.
 */
class SpinResolvedOperatorString {
private:
    // The vector representing the indices of the annihilation or creation operators in the operator string. Each index is paired with a certain spin.
    std::vector<std::pair<size_t, GQCP::Spin>> index_spin_pairs;

public:
    /*
     * MARK: Constructors
     */

    /**
     *  Construct a `SpinResolvedOperatorString` from the vector of indices that it encapsulates, together with a vector containing the spins of the individual operators.
     * 
     *  @param index_vector                The vector containing the alpha operator indices.
     *  @param spin_vector                 The vector containing the spins associated with each individual operator index.
     * 
     * Note: The index vector element zero will be matched with the spin vector element zero, the index vector element one will be matched with spin vector element one and so on.
     */
    SpinResolvedOperatorString(const std::vector<size_t>& index_vector, const std::vector<GQCP::Spin>& spin_vector) {

        // Throw an exceptioon if the vector dimensions don't match.
        if (index_vector.size() != spin_vector.size()) {
            throw std::invalid_argument("SpinUnresolvedOperatorString(const std::vector<size_t>& index_vector, const std::vector<GQCP::Spin>& spin_vector): The dimensions of the argument vectors don't match. They should be the same.");
        }

        // Initialize the empty spin pair vector.
        auto pair_vector = std::vector<std::pair<size_t, GQCP::Spin>> {};

        // Create the pairs and fill the pair vector.
        for (int i = 0; i < index_vector.size(); i++) {
            auto pair = std::pair<size_t, GQCP::Spin> {index_vector[i], spin_vector[i]};
            pair_vector.push_back(pair);
        }

        this->index_spin_pairs = pair_vector;
    };


    /*
     *  MARK: General information
     */

    /**
     *  Retrieve the operator indices from the `SpinUnresolvedOPeratorString`.
     */
    std::vector<size_t> operatorIndices() const {
        // for each pair in the operator string, save the index and return the vector containing them.
        auto index_vector = std::vector<size_t> {};

        for (int i = 0; i < this->index_spin_pairs.size(); i++) {
            index_vector.push_back(this->index_spin_pairs[i].first);
        }

        return index_vector;
    }


    /**
     *  Retrieve the operator indices from the `SpinUnresolvedOPeratorString`.
     */
    std::vector<GQCP::Spin> operatorSpins() const {
        // for each pair in the operator string, save the index and return the vector containing them.
        auto spin_vector = std::vector<GQCP::Spin> {};

        for (int i = 0; i < this->index_spin_pairs.size(); i++) {
            spin_vector.push_back(this->index_spin_pairs[i].second);
        }

        return spin_vector;
    }


    /*
     *  MARK: Component acces
     */

    /**
     * Return the spin resolved operator string split in its two separate components, one containing the alpha operators, one containing the beta operators.
     */
    SpinResolved<SpinUnresolvedOperatorString> SpinResolve() const {
        // Split the pairs in two separate spinUnresolvedOPeratorStrings
        auto alpha_index_vector = std::vector<size_t> {};
        auto beta_index_vector = std::vector<size_t> {};

        for (int i = 0; i < this->index_spin_pairs.size(); i++) {
            if (this->index_spin_pairs[i].second == GQCP::Spin::alpha) {
                alpha_index_vector.push_back(index_spin_pairs[i].first);
            } else {
                beta_index_vector.push_back(index_spin_pairs[i].first);
            }
        }

        // Determine the phase factor needed to correspond with the fermion anti-commutation rules.
    }
};


}  // namespace GQCP
