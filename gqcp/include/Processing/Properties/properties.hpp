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


#include "ONVBasis/WaveFunction/WaveFunction.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Processing/RDM/OneRDM.hpp"


namespace GQCP {


/**
 *  @param dipole_operator          the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm                  the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Vector<double, 3> calculateElectronicDipoleMoment(const VectorSQOneElectronOperator<double>& dipole_operator, const OneRDM<double>& one_rdm);

/**
 *  Calculate the electric polarizability from the linear wave function response
 * 
 *  @param F_p          the electric response force (d^2E/dFdp)
 *  @param response     the linear wave function response
 * 
 *  @return the components of the electric polarizability
 */
Matrix<double, 3, 3> calculateElectricPolarizability(const Matrix<double, Dynamic, 3>& F_p, const Matrix<double, Dynamic, 3>& response);

/**
 *  Calculate the Dyson 'amplitudes' (the coefficients of a Dyson orbital) between two wave function expressed in the same spinor basis 
 * 
 *  @param wavefunction1        a wave function in a product ONV basis  
 *  @param wavefunction2        a wave function in a product ONV basis containing one fewer electron and the same amount of orbitals that is expressed in the same basis
 *
 *  @return a vector with the Dyson orbital amplitudes  
 */
VectorX<double> calculateDysonOrbitalCoefficients(const WaveFunction& wavefunction1, const WaveFunction& wavefunction2);


}  // namespace GQCP
