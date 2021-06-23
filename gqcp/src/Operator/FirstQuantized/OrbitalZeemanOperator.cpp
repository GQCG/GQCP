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

#include "Operator/FirstQuantized/OrbitalZeemanOperator.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */


/**
 *  Construct a `OrbitalZeemanOperator` from its underlying homogeneous magnetic field and a reference point.
 * 
 *  @param B                    The external, homogeneous magnetic field.
 *  @param reference            The reference point about which the orbital Zeeman operator is calculated.
 */
OrbitalZeemanOperator::OrbitalZeemanOperator(const HomogeneousMagneticField& B, const Vector<double, 3>& reference) :
    B {B},
    BaseReferenceDependentOperator(reference) {}


/**
 *  Construct a `OrbitalZeemanOperator` from its underlying homogeneous magnetic field.
 * 
 *  @param B                    The external, homogeneous magnetic field.
 * 
 *  @note The reference point is chosen to be the gauge origin of the magnetic field.
 */
OrbitalZeemanOperator::OrbitalZeemanOperator(const HomogeneousMagneticField& B) :
    OrbitalZeemanOperator(B, B.gaugeOrigin()) {}


}  // namespace GQCP
