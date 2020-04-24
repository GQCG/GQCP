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

#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreDOCI.hpp"

#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param onv_basis       the frozen ONV basis (identical for alpha and beta)
 */
FrozenCoreDOCI::FrozenCoreDOCI(const SpinUnresolvedFrozenONVBasis& onv_basis) :
    FrozenCoreCI(std::make_shared<DOCI>(onv_basis.get_active_fock_space()), onv_basis.get_number_of_frozen_orbitals()),
    onv_basis {onv_basis} {}


}  // namespace GQCP
