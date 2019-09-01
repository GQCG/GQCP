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
 *  @return a KineticOperator
 */
KineticOperator Operator::Kinetic() {

    return KineticOperator();
}


/**
 *  @param mol              the molecule that contains the nuclear framework
 * 
 *  @return a NuclearAttractionOperator
 */
NuclearAttractionOperator Operator::NuclearAttraction(const Molecule& mol) {

    return NuclearAttractionOperator(mol.nuclearFramework());
}


/**
 *  @param origin               the origin of the dipole operator
 * 
 *  @return an ElectronicDipoleOperator
 */
ElectronicDipoleOperator Operator::ElectronicDipole(const Vector<double, 3>& o) {

    return ElectronicDipoleOperator(o);
}


/**
 *  @return a CoulombRepulsionOperator
 */
CoulombRepulsionOperator Operator::Coulomb() {

    return CoulombRepulsionOperator();
}


/**
 *  @param mol              the molecule that contains the nuclear framework
 * 
 *  @return a NuclearRepulsionOperator
 */
NuclearRepulsionOperator Operator::NuclearRepulsion(const Molecule& mol) {

    return NuclearRepulsionOperator(mol.nuclearFramework());
}


/**
 *  @param mol                  the molecule that contains the nuclear framework
 *  @param origin               the origin of the dipole operator
 * 
 *  @return a NuclearD
 */
NuclearDipoleOperator Operator::NuclearDipole(const Molecule& mol, const Vector<double, 3>& o) {

    return NuclearDipoleOperator(mol.nuclearFramework(), o);
}



}  // namespace GQCP
