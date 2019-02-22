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
#ifndef RMP2_hpp
#define RMP2_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "RHF/RHF.hpp"


namespace GQCP {


/**
 *  @param ham_par      Hamiltonian parameters in an orthornomal orbital basis
 *  @param molecule     the molecule for which the energy correction should be calculated
 *  @param rhf          the converged solution to the RHF SCF equations
 *
 *  @return the RMP2 energy correction
 */
double calculateRMP2EnergyCorrection(const HamiltonianParameters<double>& ham_par, const Molecule& molecule, const RHF& rhf);


}  // namespace GQCP



#endif /* RMP2_hpp */
