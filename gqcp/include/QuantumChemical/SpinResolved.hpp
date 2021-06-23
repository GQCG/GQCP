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


#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  A utility type encapsulating an alpha- and beta-type. In contrast to `SpinResolvedBase`, this class is used as an instantiatable type, much like you would use a std::vector<_Of>.
 * 
 *  @param _Of          The type of the alpha- and beta-objects. The name 'Of' is chosen for a natural reading `SpinResolved<_Of>`.
 */
template <typename _Of>
class SpinResolved:
    public SpinResolvedBase<_Of, SpinResolved<_Of>> {
public:
    // The type of the alpha- and beta-objects. The name 'Of' is chosen for a natural reading `SpinResolved<_Of, _Derived>`.
    using Of = _Of;
    using ComponentType = Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<Of, SpinResolved<Of>>::SpinResolvedBase;
};


}  // namespace GQCP
