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


#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/RHF/RHF.hpp"


namespace GQCP {


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthornomal basis
 *  @param molecule             the molecule for which the energy correction should be calculated
 *  @param rhf                  the converged solution to the RHF SCF equations
 *
 *  @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const SQHamiltonian<double>& sq_hamiltonian, const Molecule& molecule, const RHF& rhf);


}  // namespace GQCP
