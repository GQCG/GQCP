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
#include "FockSpace/WaveFunction/WaveFunction.hpp"
#include "Mathematical/Representation/Matrix.hpp"



namespace GQCP {


/**
 *  A class that represents a selected wave function: (potentially non normalized) expansion coefficients in a Fock space
 */
class SelectedWaveFunction {
private:
    SelectedFockSpace fock_space;
    VectorX<double> coefficients;  // the expansion coefficients for the selection in a Fock space

public:
    // CONSTRUCTORS
    SelectedWaveFunction() = default;

    /**
     *  Construct a selected wave function
     *
     *  @param fock_space           the selected Fock space in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    SelectedWaveFunction(const SelectedFockSpace& fock_space, const VectorX<double>& coefficients);


    /**
     *  Construct a selected wave function by selecting the configurations with the highest contribution from a different wave function
     *
     *  @param wave_function                           the wave_function to be selected from
     *  @param number_of_largest_contributions         the amount of selections made
     */
    SelectedWaveFunction(const WaveFunction& wave_function, size_t number_of_largest_contributions);


    // PUBLIC METHODS
    /**
     *  Print this wave function to an output in the GAMESS-US Expansion format
     *
     *  @param output_stream        the stream used for outputting
     */
    void printGamessExpansion(std::ostream& output_stream = std::cout) const;


    // GETTERS
    const VectorX<double>& get_coefficients() const { return coefficients; }
    const SelectedFockSpace& fockSpace() const { return fock_space; }


};


}  // namespace GQCP
