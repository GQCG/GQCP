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


#include "Basis/Transformations/JacobiRotation.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  A class that encapsulates an alpha- and a beta- Jacobirotation.
 */
class UJacobiRotation:
    public SpinResolvedBase<JacobiRotation, UJacobiRotation> {
public:
    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<JacobiRotation, UJacobiRotation>::Of;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<JacobiRotation, UJacobiRotation>::SpinResolvedBase;
};


}  // namespace GQCP
