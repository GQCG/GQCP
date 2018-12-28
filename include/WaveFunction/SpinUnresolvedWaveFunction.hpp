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
#ifndef GQCP_SPINUNRESOLVEDWAVEFUNCTION_HPP
#define GQCP_SPINUNRESOLVEDWAVEFUNCTION_HPP


#include "FockSpace/FockSpace.hpp"
#include "WaveFunction/WaveFunction.hpp"


namespace GQCP {


/**
 *  A class that represents a wave function: expansion coefficients in a (spin-orbital) Fock space
 */
class SpinUnresolvedWaveFunction : public WaveFunction {
    /**
     *  @param fock_space           the Fock space in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    SpinUnresolvedWaveFunction(const FockSpace& fock_space, const Eigen::VectorXd& coefficients);
};


}  // namespace GQCP


#endif  // GQCP_SPINUNRESOLVEDWAVEFUNCTION_HPP
