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
 *  NUCLEAR REPULSION OPERATOR - PUBLIC METHODS
 */

/**
 *  @return the scalar value of this nuclear repulsion operator
 */
double NuclearRepulsionOperator::value() const {

    double value {0.0};

    // // Sum over every unique nucleus pair
    // auto natoms = this->numberOfAtoms();
    // for (size_t i = 0; i < natoms; i++) {
    //     for (size_t j = i + 1; j < natoms; j++ ) {
    //         const auto atom1 = this->atoms[i];
    //         const auto atom2 = this->atoms[j];

    //         // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
    //         value += atom1.atomic_number * atom2.atomic_number / this->calculateInternuclearDistance(i, j);
    //     }
    // }

    return value;
}



/**
 *  @return the value of this nuclear dipole operator
 */
Vector<double, 3> NuclearRepulsionOperator::value() const {

    
}






/*
 *  OPERATOR - PUBLIC STATIC METHODS
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

    return NuclearAttractionOperator(mol.get_atoms());
}


/**
 *  @return a CoulombRepulsionOperator
 */
CoulombRepulsionOperator Operator::Coulomb() {

    return CoulombRepulsionOperator();
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
