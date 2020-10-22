// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDensityOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/LinearMomentumOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"


namespace GQCP {


/**
 *  A factory class for easily creating various types of first-quantized operators.
 */
class Operator {
public:
    /*
     * MARK: Scalar one-electron operators
     */

    /**
     *  Create an ElectronicDensityOperator.
     * 
     *  @return An ElectronicDensityOperator.
     */
    static ElectronicDensityOperator ElectronicDensity() { return ElectronicDensityOperator(); }

    /**
     *  Create a KineticOperator.
     * 
     *  @return A KineticOperator.
     */
    static KineticOperator Kinetic() { return KineticOperator(); }

    /**
     *  Create a NuclearAttractionOperator related to a nuclear framework.
     * 
     *  @param nuclear_framework            The nuclear framework.
     * 
     *  @return A NuclearAttractionOperator.
     */
    static NuclearAttractionOperator NuclearAttraction(const NuclearFramework& nuclear_framework) { return NuclearAttractionOperator(nuclear_framework); }

    /**
     *  Create a NuclearAttractionOperator related to a molecule.
     * 
     *  @param molecule              The molecule that contains the nuclear framework.
     * 
     *  @return A NuclearAttractionOperator.
     */
    static NuclearAttractionOperator NuclearAttraction(const Molecule& molecule) { return Operator::NuclearAttraction(molecule.nuclearFramework()); }

    /**
     *  Create an OverlapOperator.
     * 
     *  @return An OverlapOperator.
     */
    static OverlapOperator Overlap() { return OverlapOperator(); }


    /*
     * MARK: Vector one-electron operators
     */

    /**
     *  Create an AngularMomentumOperator about a reference point.
     * 
     *  @param reference                The reference point about which the angular momentum is defined.
     * 
     *  @return An AngularMomentumOperator.
     */
    static AngularMomentumOperator AngularMomentum(const Vector<double, 3>& reference = Vector<double, 3>::Zero()) { return AngularMomentumOperator(reference); }

    /**
     *  Create an ElectronicDipoleOperator with the given point as a reference.
     * 
     *  @param origin               The origin of the dipole operator.
     * 
     *  @return An ElectronicDipoleOperator.
     */
    static ElectronicDipoleOperator ElectronicDipole(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return ElectronicDipoleOperator(origin); }

    /**
     *  Create an ElectronicSpinOperator.
     * 
     *  @return An ElectronicSpinOperator.
     */
    static ElectronicSpinOperator ElectronicSpin() { return ElectronicSpinOperator(); }

    /**
     *  Create a LinearMomentumOperator.
     * 
     *  @return a LinearMomentumOperator.
     */
    static LinearMomentumOperator LinearMomentum() { return LinearMomentumOperator(); }


    /*
     *  MARK: Two-electron operators
     */

    /**
     *  Create a CoulombRepulsionOperator.
     * 
     *  @return A CoulombRepulsionOperator.
     */
    static CoulombRepulsionOperator Coulomb() { return CoulombRepulsionOperator(); }


    /*
     * MARK: Nuclear operators
    */

    /**
     *  Create a NuclearDipoleOperator from a molecule and a reference point for the dipole.
     * 
     *  @param molecule             The molecule that contains the nuclear framework.
     *  @param origin               The origin of the dipole operator.
     * 
     *  @return A NuclearDipoleOperator.
     */
    static NuclearDipoleOperator NuclearDipole(const Molecule& molecule, const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return NuclearDipoleOperator(molecule.nuclearFramework(), origin); }

    /**
     *  Create a NuclearRepulsionOperator from a molecule.
     * 
     *  @param molecule              The molecule that contains the nuclear framework.
     * 
     *  @return A NuclearRepulsionOperator.
     */
    static NuclearRepulsionOperator NuclearRepulsion(const Molecule& molecule) { return NuclearRepulsionOperator(molecule.nuclearFramework()); }
};


}  // namespace GQCP
