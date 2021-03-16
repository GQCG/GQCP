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
#include "Operator/FirstQuantized/ElectronicQuadrupoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpin_zOperator.hpp"
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
     *  @return An `ElectronicDensityOperator`.
     */
    static ElectronicDensityOperator ElectronicDensity() { return ElectronicDensityOperator(); }

    /**
     *  @return A `KineticOperator`.
     */
    static KineticOperator Kinetic() { return KineticOperator(); }

    /**
     *  @param nuclear_framework            The nuclear framework.
     * 
     *  @return A `NuclearAttractionOperator` corresponding to the given nuclear framework.
     */
    static NuclearAttractionOperator NuclearAttraction(const NuclearFramework& nuclear_framework) { return NuclearAttractionOperator(nuclear_framework); }

    /**
     *  @param molecule              The molecule that contains the nuclear framework.
     * 
     *  @return A `NuclearAttractionOperator` corresponding to the given molecule.
     */
    static NuclearAttractionOperator NuclearAttraction(const Molecule& molecule) { return Operator::NuclearAttraction(molecule.nuclearFramework()); }

    /**
     *  @return An `OverlapOperator`.
     */
    static OverlapOperator Overlap() { return OverlapOperator(); }

    /**
     *  @return An `ElectronicSpin_zOperator`.
     */
    static ElectronicSpin_zOperator ElectronicSpin_z() { return ElectronicSpin_zOperator(); }


    /*
     * MARK: Vector one-electron operators
     */

    /**
     *  @param reference                The reference point about which the angular momentum is defined.
     * 
     *  @return An `AngularMomentumOperator` relative to the given reference point.
     */
    static AngularMomentumOperator AngularMomentum(const Vector<double, 3>& reference = Vector<double, 3>::Zero()) { return AngularMomentumOperator(reference); }

    /**
     *  @param origin               The origin of the dipole operator.
     * 
     *  @return An `ElectronicDipoleOperator` relative to the given origin.
     */
    static ElectronicDipoleOperator ElectronicDipole(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return ElectronicDipoleOperator(origin); }

    /**
     *  @return An `ElectronicSpinOperator`.
     */
    static ElectronicSpinOperator ElectronicSpin() { return ElectronicSpinOperator(); }

    /**
     *  @return a `LinearMomentumOperator`.
     */
    static LinearMomentumOperator LinearMomentum() { return LinearMomentumOperator(); }


    /*
     *  MARK: Matrix one-electron operators
     */

    /**
     *  @param origin               The origin of the quadrupole operator.
     * 
     *  @return An `ElectronicQuadrupoleOperator` relative to the given origin.
     */
    static ElectronicQuadrupoleOperator ElectronicQuadrupole(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return ElectronicQuadrupoleOperator(origin); }


    /*
     *  MARK: Two-electron operators
     */

    /**
     *  @return A `CoulombRepulsionOperator`.
     */
    static CoulombRepulsionOperator Coulomb() { return CoulombRepulsionOperator(); }


    /*
     * MARK: Nuclear operators
     */

    /**
     *  @param molecule             The molecule that contains the nuclear framework.
     *  @param origin               The origin of the dipole operator.
     * 
     *  @return A `NuclearDipoleOperator` for the given molecule and reference point.
     */
    static NuclearDipoleOperator NuclearDipole(const Molecule& molecule, const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return NuclearDipoleOperator(molecule.nuclearFramework(), origin); }

    /**
     *  @param molecule              The molecule that contains the nuclear framework.
     * 
     *  @return A `NuclearRepulsionOperator` related to the given molecule.
     */
    static NuclearRepulsionOperator NuclearRepulsion(const Molecule& molecule) { return NuclearRepulsionOperator(molecule.nuclearFramework()); }
};


}  // namespace GQCP
