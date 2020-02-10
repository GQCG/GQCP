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
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreFCI.hpp"

#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param onv_basis       the spin-resolved frozen ONV basis
 */
FrozenCoreFCI::FrozenCoreFCI(const SpinResolvedFrozenONVBasis& onv_basis) :
    FrozenCoreCI(std::make_shared<FCI>(onv_basis.get_active_product_fock_space()), onv_basis.get_number_of_frozen_orbitals()),
    onv_basis (onv_basis)
{}



}  // namespace GQCP
