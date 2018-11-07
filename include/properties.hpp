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
#ifndef properties_hpp
#define properties_hpp

#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  Calculate the electronic energy as a result of the contraction of the 1- and 2-RDMs with the one- and two-electron integrals
 *
 *  @param ham_par      the Hamiltonian parameters containing the one- and two-electron integrals
 *  @param one_rdm      the 1-RDM
 *  @param two_rdm      the 2-RDM
 *
 *  @return the electronic energy
 */
double calculateElectronicEnergy(const GQCP::HamiltonianParameters& ham_par, const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm);


/**
 *  Calculate the electronic electric dipole moment
 *
 *  @param dipole_operator      the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm              the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Eigen::Vector3d calculateElectronicDipoleMoment(const std::array<GQCP::OneElectronOperator, 3>& dipole_operator, const GQCP::OneRDM& one_rdm);


}  // namespace GQCP


#endif /* properties_hpp */
