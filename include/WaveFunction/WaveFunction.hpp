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
#ifndef GQCP_WAVEFUNCTION_HPP
#define GQCP_WAVEFUNCTION_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "Mathematical/Matrix.hpp"



namespace GQCP {


/**
 *  A class that represents a wave function: expansion coefficients in a Fock space
 */
class WaveFunction {
private:
    std::shared_ptr<BaseFockSpace> fock_space;
    VectorX<double> coefficients;  // the expansion coefficients in the Fock space

public:
    // CONSTRUCTORS
    WaveFunction() = default;

    /**
     *  Construct a normalized wave function from possibly non-normalized coefficients
     *
     *  @param base_fock_space      the Fock space in which the wave function 'lives'
     *  @param coefficients         the expansion coefficients
     */
    WaveFunction(const BaseFockSpace& base_fock_space, const VectorX<double>& coefficients);


    // GETTERS
    const VectorX<double>& get_coefficients() const { return coefficients; }
    const BaseFockSpace& get_fock_space() const { return *fock_space; }


    // PUBLIC METHODS
    /**
     *  @return the Shannon entropy (or information content) of the wave function
     */
    double calculateShannonEntropy() const;

    /**
     *  Transform the underlying ONV basis of the wave function (only for FCI [ProductFockSpace])
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
     void basisTransform(const SquareMatrix<double>& T);
};


}  // namespace GQCP


#endif  // GQCP_WAVEFUNCTION_HPP
