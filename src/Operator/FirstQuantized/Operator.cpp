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

    const auto& nuclei = this->nuclear_framework.nucleiAsVector();

    // Sum over every unique nucleus pair
    double value {0.0};
    const auto n_nuclei = this->nuclearFramework().numberOfNuclei();
    for (size_t i = 0; i < n_nuclei; i++) {
        for (size_t j = i + 1; j < n_nuclei; j++ ) {
            const auto nucleus1 = nuclei[i];
            const auto nucleus2 = nuclei[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
            value += nucleus1.charge() * nucleus2.charge() / nucleus1.calculateDistance(nucleus2);
        }
    }

    return value;
}



/*
 *  NUCLEAR DIPOLE OPERATOR - CONSTRUCTOR
 */

/**
 *  @param nuclear_framework            the nuclear framework underlying a nuclear operator
 *  @param o                            the origin of the multipole
 */
NuclearDipoleOperator::NuclearDipoleOperator(const NuclearFramework& nuclear_framework, const Vector<double, 3>& o) :
    BaseNuclearOperator(nuclear_framework),
    BaseMultipoleOperator(o)
{}



/*
 *  NUCLEAR DIPOLE OPERATOR - PUBLIC METHODS
 */

/**
 *  @return the value of this nuclear dipole operator
 */
Vector<double, 3> NuclearDipoleOperator::value() const {

    Vector<double, 3> mu = Vector<double, 3>::Zero();

    const auto& nuclei = this->nuclear_framework.nucleiAsVector();
    for (const auto& nucleus : nuclei) {
        mu += nucleus.charge() * nucleus.position();
    }

    return mu;
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
