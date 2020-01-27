// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#include "FockSpace/WaveFunction/WaveFunctionSelection.hpp"

#include <queue>


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Construct a selected wave function
 *
 *  @param fock_space           the selected Fock space in which the wave function 'lives'
 *  @param coefficients         the expansion coefficients
 */
WaveFunctionSelection::WaveFunctionSelection(const SelectedFockSpace& fock_space, const VectorX<double>& coefficients) :
    configurations (fock_space.get_configurations()),
    coefficients (coefficients)
{
}



/**
 *  Construct a selected wave function by selecting the configurations with the highest contribution from a different wave function
 *
 *  @param wave_function                           the wave function to be selected from
 *  @param number_of_highest_contributions         the amount of selections made
 */
WaveFunctionSelection::WaveFunctionSelection(const WaveFunction& wave_function, size_t number_of_largest_contributions) {

    // A struct containing a coefficient with its address in a Fock space
    struct AddressCoefficientPair 
    {
        double coeff;
        size_t address;
    };

    // Struct for the priority queue API to compare the coefficients of two AddressCoefficientPairs
    struct AddressCoefficientPairComparer 
    { 
        int operator() (const AddressCoefficientPair& p1, const AddressCoefficientPair& p2) 
        { 
            return std::abs(p1.coeff) >std::abs(p2.coeff); 
        } 
    }; 

    // Initialize a "min heap" datastructure for which the top value is the lowest value
    std::priority_queue <AddressCoefficientPair, std::vector<AddressCoefficientPair>, AddressCoefficientPairComparer > min_heap; 

    // Extract the amount of requested largest contributions
    const auto coefficients = wave_function.get_coefficients();
    for (size_t i = 0; i < wave_function.get_fock_space().get_dimension(); i++) {
        min_heap.push(AddressCoefficientPair {coefficients(i), i});
        // When we've reached the requested amount of contributions start to remove the top element from the heap (the one with the smallest contribution)
        if (min_heap.size() > number_of_largest_contributions) {
            min_heap.pop();
        }
    }

    this->coefficients = VectorX<double>(number_of_largest_contributions);
    std::vector<Configuration> configurations;
    configurations.reserve(number_of_largest_contributions);

    // Generate the selection from the stored addresses
    const BaseFockSpace& fock_space = wave_function.get_fock_space();

    size_t i = 0;
    while (!min_heap.empty()) {
        const AddressCoefficientPair& x = min_heap.top();
        this->coefficients (i) = x.coeff;
        configurations.push_back(fock_space.configuration(x.address));
        min_heap.pop();
        i++;
    }   
    
    this->configurations = configurations;
}


// PUBLIC METHODS
/**
 *  Print this wave function to an output in the GAMESS-US Expansion format
 *
 *  @param output_stream        the stream used for outputting
 */
void WaveFunctionSelection::printGamessExpansion(std::ostream& output_stream) const {
    for (size_t i = 0; i < configurations.size(); i++) {
        output_stream << configurations[i].asCompactString(" ") << " | " << this->coefficients(i) << "\n";
    }
}


}  // namespace GQCP
