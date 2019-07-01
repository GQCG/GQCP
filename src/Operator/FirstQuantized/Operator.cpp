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
#include "Operator/FirstQuantized/Operator.hpp"


namespace GQCP {


/*
 *  PUBLIC STATIC METHODS
 */

/**
 *  @return an OverlapOperator
 */
OverlapOperator Operator::Overlap() {

    return OverlapOperator();
}


/**
 *  @return a KineticEnergyOperator
 */
KineticEnergyOperator Operator::Kinetic() {

    return KineticEnergyOperator();
}


/**
 *  @param mol              the molecule that contains the nuclear framework
 * 
 *  @return a NuclearAttractionEnergyOperator
 */
NuclearAttractionEnergyOperator Operator::NuclearAttraction(const Molecule& mol) {

    return NuclearAttractionEnergyOperator(mol.get_atoms());
}


/**
 *  @return a CoulombInteractionEnergyOperator
 */
CoulombInteractionEnergyOperator Operator::Coulomb() {

    return CoulombInteractionEnergyOperator();
}


/**
 *  @param origin               the origin of the dipole operator
 * 
 *  @return an ElectronicDipoleOperator
 */
ElectronicDipoleOperator Operator::ElectronicDipole(const Vector<double, 3>& o = Vector<double, 3>::Zero(3)) {

    return ElectronicDipoleOperator(o);
}


/**
 *  @param mol                  the molecule that contains the nuclear framework
 *  @param origin               the origin of the dipole operator
 * 
 *  @return a NuclearD
 */
NuclearDipoleOperator Operator::NuclearDipole(const Molecule& mol, const Vector<double, 3>& o = Vector<double, 3>::Zero(3)) {

    return NuclearDipoleOperator(mol.get_atoms(), o);
}



}  // namespace GQCP
