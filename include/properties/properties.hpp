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
#ifndef properties_hpp
#define properties_hpp

#include "Operator/OneElectronOperator.hpp"
#include "RDM/OneRDM.hpp"


namespace GQCP {


/**
 *  @param dipole_operator      the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm              the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Eigen::Vector3d calculateElectronicDipoleMoment(const std::array<OneElectronOperator<double>, 3>& dipole_operator, const OneRDM& one_rdm);



}  // namespace GQCP


#endif /* properties_hpp */
