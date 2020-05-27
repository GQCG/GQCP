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

#include "Processing/RDM/FrozenCoreFCIRDMBuilder.hpp"

#include "Processing/RDM/FCIRDMBuilder.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param onv_basis       the frozen spin-resolved ONV basis
 */
FrozenCoreFCIRDMBuilder::FrozenCoreFCIRDMBuilder(const SpinResolvedFrozenONVBasis& onv_basis) :
    FrozenCoreRDMBuilder(std::make_shared<FCIRDMBuilder>(onv_basis.activeONVBasis()), onv_basis.numberOfFrozenOrbitals()),
    onv_basis {onv_basis} {}


}  // namespace GQCP
