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
#include "HamiltonianBuilder/FrozenCoreDOCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the frozen Fock space (identical for alpha and beta)
 */
FrozenCoreDOCI::FrozenCoreDOCI(const FrozenFockSpace& fock_space) :
    FrozenCoreCI(std::make_shared<DOCI>(fock_space.get_active_fock_space()), fock_space.get_number_of_frozen_orbitals()),
    fock_space (fock_space)
{}



}  // namespace GQCP
