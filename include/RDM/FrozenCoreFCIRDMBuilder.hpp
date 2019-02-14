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
#ifndef GQCP_FROZENCOREFCIRDMBUILDER_HPP
#define GQCP_FROZENCOREFCIRDMBUILDER_HPP


#include "FockSpace/FrozenProductFockSpace.hpp"
#include "RDM/FrozenCoreRDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-RDMs from wave functions expanded in the frozen full CI product Fock space
 */
class FrozenCoreFCIRDMBuilder : public FrozenCoreRDMBuilder {
    FrozenProductFockSpace fock_space;  // contains both the frozen alpha and beta Fock space

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the frozen product Fock space
     */
    explicit FrozenCoreFCIRDMBuilder(const FrozenProductFockSpace& fock_space);

    // OVERRIDDEN GETTERS
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_FROZENCOREFCIRDMBUILDER_HPP
