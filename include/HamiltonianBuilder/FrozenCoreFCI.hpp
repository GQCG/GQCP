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
#ifndef GQCP_FROZENCOREFCI_HPP
#define GQCP_FROZENCOREFCI_HPP


#include "HamiltonianBuilder/FrozenCoreCI.hpp"
#include "FockSpace/FrozenProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"


namespace GQCP {


/**
 *  A class capable of generating the matrix representation of the frozen core FCI Hamiltonian
 */
class FrozenCoreFCI : public FrozenCoreCI {
private:
    FrozenProductFockSpace fock_space;  // contains both the frozen alpha and beta Fock space

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the frozen product Fock space
     */
    explicit FrozenCoreFCI(const FrozenProductFockSpace& fock_space);


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_FCI_HPP
