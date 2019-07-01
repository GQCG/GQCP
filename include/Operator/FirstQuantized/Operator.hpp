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
#ifndef GQCP_OPERATOR_HPP
#define GQCP_OPERATOR_HPP


#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"
#include "Operator/FirstQuantized/BaseMultipoleOperator.hpp"
#include "Molecule.hpp"



namespace GQCP {


/*
 *  OPERATORS 
 */


/**
 *  A class that represents the overlap operator
 */
class OverlapOperator {};


/**
 *  A class that represents the kinetic energy operator for the electrons
 */
class KineticEnergyOperator {};


/**
 *  A class that represents the nuclear attraction energy operator for the electrons
 */
class NuclearAttractionEnergyOperator : public BaseNuclearOperator {
public:
    using BaseNuclearOperator::BaseNuclearOperator;  // inherit base constructors
};


/**
 *  A class that represents the Coulomb interaction energy operator between the electrons
 */
class CoulombInteractionEnergyOperator {};


/**
 *  A class that represents the electronic dipole operator for the electrons
 */
class ElectronicDipoleOperator: public BaseMultipoleOperator {
public:
    using BaseMultipoleOperator::BaseMultipoleOperator;  // inherit base constructors
};


/**
 *  A class that represents the nuclear dipole operator
 */
class NuclearDipoleOperator: public BaseNuclearOperator, public BaseMultipoleOperator {
public:
    // CONSTRUCTORS

    /**
     *  @param atoms            the atoms that represent the nuclear framework
     *  @param o                the origin of the multipole
     */
    NuclearDipoleOperator(const std::vector<Atom>& atoms, const Vector<double, 3>& o = Vector<double, 3>::Zero(3));
};



/*
 *  OPERATOR
 */

/**
 *  A class that is used to construct operators, much like a factory class does
 */
class Operator {
public:

    // PUBLIC STATIC METHODS

    /**
     *  @return an OverlapOperator
     */
    static OverlapOperator Overlap();

    /**
     *  @return a KineticEnergyOperator
     */
    static KineticEnergyOperator Kinetic();

    /**
     *  @param mol              the molecule that contains the nuclear framework
     * 
     *  @return a NuclearAttractionEnergyOperator
     */
    static NuclearAttractionEnergyOperator NuclearAttraction(const Molecule& mol);

    /**
     *  @return a CoulombInteractionEnergyOperator
     */
    static CoulombInteractionEnergyOperator Coulomb();

    /**
     *  @param origin               the origin of the dipole operator
     * 
     *  @return an ElectronicDipoleOperator
     */
    static ElectronicDipoleOperator ElectronicDipole(const Vector<double, 3>& o = Vector<double, 3>::Zero(3));

    /**
     *  @param mol                  the molecule that contains the nuclear framework
     *  @param origin               the origin of the dipole operator
     * 
     *  @return a NuclearDipoleOperator
     */
    static NuclearDipoleOperator NuclearDipole(const Molecule& mol, const Vector<double, 3>& o = Vector<double, 3>::Zero(3));
};


}  // namespace GQCP



#endif  // GQCP_OPERATOR_HPP
