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
#include "Mathematical/Representation/BlockRankFourTensor.hpp"
#include "Mathematical/Representation/Tensor.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/OneRDM.hpp"


namespace GQCP {


/**
 *  A class that represents a converged solution to the RHF SCF equations
 */
class RHF {
private:
    double electronic_energy;
    TransformationMatrix<double> C;  // transformation matrix from the AO basis to the RHF MO basis
    VectorX<double> orbital_energies;  // sorted in ascending energies


public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    RHF();  // need default constructor

    /**
     *  Constructor based on given converged solutions of the RHF SCF equations
     *
     *  @param electronic_energy    the converged RHF electronic energy
     *  @param C                    the coefficient matrix, i.e. the transformation matrix from the AO basis to the RHF MO basis
     *  @param orbital_energies     the RHF MO energies
     */
    RHF(double electronic_energy, const TransformationMatrix<double>& C, const VectorX<double>& orbital_energies);


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    const TransformationMatrix<double>& get_C() const { return this->C; }
    const VectorX<double>& get_orbital_energies() const { return this->orbital_energies; }
    double get_orbital_energies(size_t index) const { return this->orbital_energies(index); }
};


}  // namespace GQCP
