// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_WAVEFUNCTION_HPP
#define GQCP_WAVEFUNCTION_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include <memory.h>
#include "common.hpp"



namespace GQCP {


/**
 *  A class that represents a wave function: expansion coefficients in a Fock space
 */
class WaveFunction {
private:
    std::shared_ptr<BaseFockSpace> fock_space;
    Eigen::VectorXd coefficients;  // Expansion coefficients of a wave function in the Fock space

public:
    // CONSTRUCTORS
    WaveFunction() = default;

    /**
     *  @param base_fock_space      the Fock space in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    WaveFunction(const BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients);


    // GETTERS
    const Eigen::VectorXd& get_coefficients() const { return coefficients; }
    const BaseFockSpace& get_fock_space() const { return *fock_space; }


    // PUBLIC METHODS
    /**
     *  @return the Shannon entropy (or information content) of the wave function
     */
    double calculateShannonEntropy() const;
};


}  // namespace GQCP


#endif  // GQCP_WAVEFUNCTION_HPP
