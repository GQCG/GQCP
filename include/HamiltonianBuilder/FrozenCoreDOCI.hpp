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
#ifndef GQCP_FROZENCOREDOCI_HPP
#define GQCP_FROZENCOREDOCI_HPP


#include "HamiltonianBuilder/FrozenCore.hpp"
#include "HamiltonianBuilder/FrozenFockSpace.hpp"
#include "HamiltonianBuilder/DOCI.hpp"


namespace GQCP {


/**
 *  A HamiltonianBuilder for frozen core DOCI: it builds the matrix representation of the frozen core DOCI Hamiltonian, in a Fock space where orbitals are either doubly occupied or unoccupied.
 */
class FrozenCoreDOCI : public FrozenCore {
private:
    FrozenFockSpace fock_space;  // both the alpha and beta Fock space

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the frozen Fock space, identical for alpha and beta
     */
    explicit FrozenCoreDOCI(const FrozenFockSpace& fock_space);

    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_DOCI_HPP
