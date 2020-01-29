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
#pragma once


#include "Basis/TransformationMatrix.hpp"
#include "FockSpace/SelectedFockSpace.hpp"
#include "Mathematical/Representation/Matrix.hpp"



namespace GQCP {


/**
 *  A class that represents a selection of a wave function
 */
class WaveFunctionSelection {
private:
    std::vector<Configuration> configurations;
    VectorX<double> coefficients;  // the coefficients for the ONV configurations

public:
    // CONSTRUCTORS
    WaveFunctionSelection() = default;

    /**
     *  Construct a selected wave function
     *
     *  @param configurations           a set of ONV configurations 
     *  @param coefficients             the expansion coefficients
     */
    WaveFunctionSelection(const std::vector<Configuration>& configurations, const VectorX<double>& coefficients);


    // PUBLIC METHODS
    /**
     *  Print this wave function to an output in the GAMESS-US Expansion format
     *
     *  @param output_stream        the stream used for outputting
     */
    void printGAMESSUSExpansion(std::ostream& output_stream = std::cout) const;


    // GETTERS
    const VectorX<double>& get_coefficients() const { return coefficients; }
    const std::vector<Configuration>& get_configurations() const { return configurations; }
};


}  // namespace GQCP
