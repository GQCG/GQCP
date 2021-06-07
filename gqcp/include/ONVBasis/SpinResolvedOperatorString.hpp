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

#include <math.h>


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
     *  @param index_vector                The vector containing the operator indices.
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
     *  Retrieve the operator indices from the `SpinUnresolvedOperatorString`.
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
     *  Retrieve the operator spins from the `SpinUnresolvedOperatorString`.
     */
    std::vector<GQCP::Spin> operatorSpins() const {
        // For each pair in the operator string, save the spin and return the vector containing them.
        auto spin_vector = std::vector<GQCP::Spin> {};

        for (int i = 0; i < this->index_spin_pairs.size(); i++) {
            spin_vector.push_back(this->index_spin_pairs[i].second);
        }

        return spin_vector;
    }


    /*
     *  MARK: Spin component acces
     */

    /**
     * Return the spin resolved operator string split in its two separate components, one containing the alpha operators, one containing the beta operators.
     */
    SpinResolved<SpinUnresolvedOperatorString> SpinResolve() const {
        // Split pairs in two separate spinUnresolvedOPeratorStrings.
        // Define a separate vector for the alpha and beta components.
        auto alpha_index_vector = std::vector<size_t> {};
        auto beta_index_vector = std::vector<size_t> {};

        // Loop over the original pair vector. If the pair contains an alpha spin, add the index to the alpha vector, if the pair contains a beta spin, add the index to the beta vector.
        for (int i = 0; i < this->index_spin_pairs.size(); i++) {
            if (this->index_spin_pairs[i].second == GQCP::Spin::alpha) {
                alpha_index_vector.push_back(index_spin_pairs[i].first);
            } else {
                beta_index_vector.push_back(index_spin_pairs[i].first);
            }
        }

        // Determine the phase factor needed to make the splitting of the original operator string in its spin components correspond with the fermion anti-commutation rules.
        // We will define an intermediate phase factor, which is updated throughout the procedure.
        int intermediate_phase_factor = 1;

        // First, we loop the spin string to check whether an element has alpha spin.
        for (int i = 0; i < this->operatorSpins().size(); i++) {

            // If it does, we continue, if it doesn't, we move on to the next pair in the vector.
            if (this->operatorSpins()[i] == GQCP::Spin::alpha) {

                // We check all elements left of our found alpha element, to see whether or not a beta element is found.
                for (int j = i; j >= 0; j--) {

                    // If we find a beta element left of the alpha element, we must swap the two elements in the vector and update the phase factor by multiplying it by -1, due to the anti-commutation rules.
                    // We don't physically swap the elements in the vector, but by simply checking the number of beta's left of the alpha in question, we mimic the swap.
                    if (this->operatorSpins()[j] == GQCP::Spin::beta) {
                        intermediate_phase_factor *= -1;
                    }
                }
            }
        }

        // The phase factor has been updated in the loops, in order to correspond to a modified operator string with all alpha operators first, followed by all beta operators.
        // The total operator string is applied on |vac_a>|vac_b>, e.g.: a_a a_a a_b a_b|vac_a>|vac_b>. But we want a_a a_a |vac_a> a_b a_b|vac_b>.
        // Now we have to apply the fermion anti-commutation rules for eacht time we move |vac_alpha> over a beta operator. This results in p^N_beta.
        int final_phase_factor = pow(intermediate_phase_factor, beta_index_vector.size());

        return SpinResolved<SpinUnresolvedOperatorString> {SpinUnresolvedOperatorString(alpha_index_vector, final_phase_factor), SpinUnresolvedOperatorString {beta_index_vector}};
    }
};


}  // namespace GQCP
