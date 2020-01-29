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


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Construct a selected wave function
 *
 *  @param configurations           a set of ONV configurations 
 *  @param coefficients             the expansion coefficients
 */
WaveFunctionSelection::WaveFunctionSelection(const std::vector<Configuration>& configurations, const VectorX<double>& coefficients) :
    configurations (configurations),
    coefficients (coefficients)
{}


// PUBLIC METHODS
/**
 *  Print this wave function to an output in the GAMESS-US Expansion format
 *
 *  @param output_stream        the stream used for outputting
 */
void WaveFunctionSelection::printGamessExpansion(std::ostream& output_stream) const {
    for (size_t i = 0; i < configurations.size(); i++) {
        output_stream << configurations[i].spinSummedRepresentation(" ") << " | " << this->coefficients(i) << "\n";
    }
}


}  // namespace GQCP
