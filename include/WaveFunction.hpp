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

#include "common.hpp"



namespace GQCP {


/**
 *  WaveFunction contains the expansion coefficients in its given FockSpace
 */
class WaveFunction {
private:
    BaseFockSpace* fock_space;
    Eigen::VectorXd coefficients;  // Expansion coefficients of a wave function in the Fock space

public:
    // CONSTRUCTORS
    WaveFunction(BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients);


    // GETTERS
    Eigen::VectorXd get_coefficients() const { return coefficients; }
    BaseFockSpace& get_fock_space() const { return *fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_WAVEFUNCTION_HPP
