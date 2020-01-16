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
#include "QCMethod/RHF/RHF.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Default constructor setting everything to zero
 */
RHF::RHF() :
    RHF(0.0, TransformationMatrix<double>::Zero(0, 0), VectorX<double>::Zero(0))
{}


/**
 *  Constructor based on given converged solutions of the RHF SCF equations
 *
 *  @param electronic_energy    the converged RHF electronic energy
 *  @param C                    the coefficient matrix, i.e. the transformation matrix from the AO basis to the RHF MO basis
 *  @param orbital_energies     the RHF MO energies
 */
RHF::RHF(double electronic_energy, const TransformationMatrix<double>& C, const VectorX<double>& orbital_energies) :
    electronic_energy (electronic_energy),
    C (C),
    orbital_energies (orbital_energies)
{}


}  // namespace GQCP
