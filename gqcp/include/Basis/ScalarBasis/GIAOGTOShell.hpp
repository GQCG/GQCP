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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Utilities/HomogeneousMagneticField.hpp"


namespace GQCP {


/**
 *  A gauge-including GTO shell.
 */
class GIAOGTOShell:
    public GTOShell {
public:
    // FIXME: For initial compilations, we didn't implement the phase factor per se for every basis function. We should implement a GIAO-modification to a CartesianGTO.
    using Primitive = GTOShell::Primitive;
    using BasisFunction = GTOShell::Primitive;


private:
    // The specificiation of the magnetic field for this GIAO.
    HomogeneousMagneticField m_field;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a gauge-including shell from a field-free shell and a magnetic field.
     * 
     *  @param shell            The field-free GTO shell.
     *  @param field            The specificiation of the magnetic field for this GIAO.
     */
    GIAOGTOShell(const GTOShell& shell, const HomogeneousMagneticField& field) :
        GTOShell(shell),
        m_field {field} {}


    /*
     *  MARK: Access
     */

    /**
     *  @return The specificiation of the magnetic field for this GIAO.
     */
    const HomogeneousMagneticField& field() const { return this->m_field; }
};


}  // namespace GQCP
