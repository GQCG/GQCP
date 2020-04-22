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
#pragma once


#include "Operator/FirstQuantized/BaseMultipoleOperator.hpp"
#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"
#include <cstddef>


namespace GQCP {


/**
 *  A class that represents the nuclear dipole operator
 */
class NuclearDipoleOperator: public BaseNuclearOperator, public BaseMultipoleOperator {
public:
    // CONSTRUCTORS

    /**
     *  @param nuclear_framework            the nuclear framework underlying a nuclear operator
     *  @param o                            the origin of the multipole
     */
    NuclearDipoleOperator(const NuclearFramework& nuclear_framework, const Vector<double, 3>& o = Vector<double, 3>::Zero());


    // PUBLIC METHODS

    /**
     *  @return the value of this nuclear dipole operator
     */
    Vector<double, 3> value() const;
};


}  // namespace GQCP
