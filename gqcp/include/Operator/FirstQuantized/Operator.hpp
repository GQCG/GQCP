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
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDensityOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
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
    // PUBLIC STATIC METHODS

    /**
     *  @return a CoulombRepulsionOperator
     */
    static CoulombRepulsionOperator Coulomb() { return CoulombRepulsionOperator(); }

    /**
     *  @return an ElectronicDensityOperator
     */
    static ElectronicDensityOperator ElectronicDensity() { return ElectronicDensityOperator(); }

    /**
     *  @param origin               the origin of the dipole operator
     * 
     *  @return an ElectronicDipoleOperator
     */
    static ElectronicDipoleOperator ElectronicDipole(const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return ElectronicDipoleOperator(origin); }

    /**
     *  @return an ElectronicSpinOperator
     */
    static ElectronicSpinOperator ElectronicSpin() { return ElectronicSpinOperator(); }

    /**
     *  @return a KineticOperator
     */
    static KineticOperator Kinetic() { return KineticOperator(); }

    /**
     *  @param nuclear_framework            the nuclear framework
     * 
     *  @return a NuclearAttractionOperator
     */
    static NuclearAttractionOperator NuclearAttraction(const NuclearFramework& nuclear_framework) { return NuclearAttractionOperator(nuclear_framework); }

    /**
     *  @param molecule              the molecule that contains the nuclear framework
     * 
     *  @return a NuclearAttractionOperator
     */
    static NuclearAttractionOperator NuclearAttraction(const Molecule& molecule) { return Operator::NuclearAttraction(molecule.nuclearFramework()); }

    /**
     *  @param molecule             the molecule that contains the nuclear framework
     *  @param origin               the origin of the dipole operator
     * 
     *  @return a NuclearDipoleOperator
     */
    static NuclearDipoleOperator NuclearDipole(const Molecule& molecule, const Vector<double, 3>& origin = Vector<double, 3>::Zero()) { return NuclearDipoleOperator(molecule.nuclearFramework(), origin); }

    /**
     *  @param molecule              the molecule that contains the nuclear framework
     * 
     *  @return a NuclearRepulsionOperator
     */
    static NuclearRepulsionOperator NuclearRepulsion(const Molecule& molecule) { return NuclearRepulsionOperator(molecule.nuclearFramework()); }

    /**
     *  @return an OverlapOperator
     */
    static OverlapOperator Overlap() { return OverlapOperator(); }
};


}  // namespace GQCP
