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
#ifndef RMP2_hpp
#define RMP2_hpp

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "RHF/RHF.hpp"


namespace GQCP {


/**
 *  @return the RMP2 energy correction based on given @param Hamiltonian parameters ham_par, a given @param molecule and a converged solution @param rhf to the RHF SCF equations
 */
double calculateRMP2EnergyCorrection(const GQCP::HamiltonianParameters& ham_par, const GQCP::Molecule& molecule, const GQCP::RHF& rhf);


}  // namespace GQCP



#endif /* RMP2_hpp */
