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

#include "ONVBasis/SpinResolvedOperatorString.hpp"


namespace GQCP {


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
SpinResolvedOperatorString::SpinResolvedOperatorString(const std::vector<size_t>& index_vector, const std::vector<Spin>& spin_vector) {

    // Throw an exception if the vector dimensions don't match.
    if (index_vector.size() != spin_vector.size()) {
        throw std::invalid_argument("SpinUnresolvedOperatorString(const std::vector<size_t>& index_vector, const std::vector<Spin>& spin_vector): The dimensions of the argument vectors don't match. They should be the same.");
    }

    // Initialize the empty spin pair vector.
    auto pair_vector = std::vector<std::pair<size_t, Spin>> {};

    // Create the pairs and fill the pair vector.
    for (int i = 0; i < index_vector.size(); i++) {
        auto pair = std::pair<size_t, Spin> {index_vector[i], spin_vector[i]};
        this->index_spin_pairs.push_back(pair);
    }
}


SpinResolvedOperatorString SpinResolvedOperatorString::FromONV(const SpinResolvedONV& onv) {

    const auto& onv_alpha = onv.onv(Spin::alpha);
    const auto& onv_beta = onv.onv(Spin::beta);

    // The occupied indices from the alpha part of a `SpinResolvedONV` can be used to represent the alpha operator indices. The same holds for the beta indices. The beta operator indices are added to the total operator string, behind the alpha indices.
    std::vector<size_t> index_vector = onv_alpha.occupiedIndices();
    index_vector.insert(index_vector.end(), onv_beta.occupiedIndices().begin(), onv_beta.occupiedIndices().end());

    // To determine the spins, we can use the number of alpha and beta electrons respectively. Once again, the beta components are appended at the end of the alpha components.
    std::vector<Spin> spin_vector {onv_alpha.numberOfElectrons(), Spin::alpha};
    spin_vector.insert(spin_vector.end(), onv_beta.numberOfElectrons(), Spin::beta);

    return SpinResolvedOperatorString(index_vector, spin_vector);
}


/*
*  MARK: General information
*/

/**
 *  Retrieve the operator indices from the `SpinResolvedOperatorString`.
 * 
 *  @return The operator indices as a vector.
 */
std::vector<size_t> SpinResolvedOperatorString::operatorIndices() const {
    // For each pair in the operator string, save the index and return the vector containing them.
    auto index_vector = std::vector<size_t> {};

    for (int i = 0; i < this->index_spin_pairs.size(); i++) {
        index_vector.push_back(this->index_spin_pairs[i].first);
    }

    return index_vector;
}


/**
 *  Retrieve the operator spins from the `SpinResolvedOperatorString`.
 * 
 *  @return The operator spins as a vector.
 */
std::vector<Spin> SpinResolvedOperatorString::operatorSpins() const {
    // For each pair in the operator string, save the spin and return the vector containing them.
    auto spin_vector = std::vector<Spin> {};

    for (int i = 0; i < this->index_spin_pairs.size(); i++) {
        spin_vector.push_back(this->index_spin_pairs[i].second);
    }

    return spin_vector;
}


/*
*  MARK: Spin component acces
*/

/**
 *  Split the spin-resolved operator string into its two separate components, one containing the alpha operators, one containing the beta operators. The alpha operator string will recieve the total phase factor assiciated with this action, the phase factor of the beta string will be one.
 * 
 *  @return The two separate components of the spin-resolved operator string.
 */
SpinResolved<SpinUnresolvedOperatorString> SpinResolvedOperatorString::spinResolve() const {
    // Split pairs in two separate spinUnresolvedOperatorStrings.
    // Define a separate vector for the alpha and beta components.
    auto alpha_index_vector = std::vector<size_t> {};
    auto beta_index_vector = std::vector<size_t> {};

    // Loop over the original pair vector. If the pair contains an alpha spin, add the index to the alpha vector, if the pair contains a beta spin, add the index to the beta vector.
    for (int i = 0; i < this->index_spin_pairs.size(); i++) {
        if (this->index_spin_pairs[i].second == Spin::alpha) {
            alpha_index_vector.push_back(index_spin_pairs[i].first);
        } else {
            beta_index_vector.push_back(index_spin_pairs[i].first);
        }
    }

    // Determine the phase factor needed to make the splitting of the original operator string in its spin components correspond with the fermion anti-commutation rules. In essence, we're counting the number of beta operators that appear in front of every alpha operator.
    // We will define an intermediate phase factor, which is updated throughout the procedure.
    int phase_factor = 1;

    // First, we loop the spin string to check whether an element has alpha spin.
    for (int i = 0; i < this->operatorSpins().size(); i++) {

        // If it does, we continue, if it doesn't, we move on to the next pair in the vector.
        if (this->operatorSpins()[i] == Spin::alpha) {

            // We check all elements left of our found alpha element, to see whether or not a beta element is found.
            for (int j = i; j >= 0; j--) {

                // If we find a beta element left of the alpha element, we must swap the two elements in the vector and update the phase factor by multiplying it by -1, due to the anti-commutation rules.
                // We don't physically swap the elements in the vector, but by simply checking the number of beta's left of the alpha in question, we mimic the swap.
                if (this->operatorSpins()[j] == Spin::beta) {
                    phase_factor *= -1;
                }
            }
        }
    }

    return SpinResolved<SpinUnresolvedOperatorString> {SpinUnresolvedOperatorString(alpha_index_vector, phase_factor), SpinUnresolvedOperatorString {beta_index_vector}};
}


}  // namespace GQCP
