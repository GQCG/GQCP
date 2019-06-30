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
     *  @param atoms                the atoms that represent the nuclear framework
     *  @param o        the origin of the multipole
     */
    NuclearDipoleOperator(const std::vector<Atom>& atoms, const Vector<double, 3>& o = Vector<double, 3>::Zero(3));
};


}  // namespace GQCP



#endif  // GQCP_OPERATOR_HPP
