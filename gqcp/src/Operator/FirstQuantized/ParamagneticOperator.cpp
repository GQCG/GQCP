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

#include "Operator/FirstQuantized/ParamagneticOperator.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param L                The angular momentum operator.
 *  @param B                The external, homogeneous magnetic field.
 */
ParamagneticOperator::ParamagneticOperator(const AngularMomentumOperator& L, const HomogeneousMagneticField& B) :
    B {B},
    L {L} {}


/**
 *  Construct a `ParamagneticOperator` from its underlying homogeneous magnetic field. The reference point for the calculation of the angular momentum operator is the gauge origin of the magnetic field.
 * 
 *  @param B                The external, homogeneous magnetic field.
 */
ParamagneticOperator::ParamagneticOperator(const HomogeneousMagneticField& B) :
    ParamagneticOperator(AngularMomentumOperator(B.gaugeOrigin()), B) {}


}  // namespace GQCP
