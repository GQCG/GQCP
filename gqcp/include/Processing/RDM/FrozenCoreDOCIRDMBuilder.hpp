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


#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "Processing/RDM/FrozenCoreRDMBuilder.hpp"
#include "Processing/RDM/RDMs.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-RDMs from wave functions expanded in the frozen DOCI ONV basis
 */
class FrozenCoreDOCIRDMBuilder: public FrozenCoreRDMBuilder {
    SpinUnresolvedFrozenONVBasis fock_space;  // both the frozen alpha and beta spin-unresolved ONV basis

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       both the frozen alpha and beta spin-unresolved ONV basis
     */
    explicit FrozenCoreDOCIRDMBuilder(const SpinUnresolvedFrozenONVBasis& fock_space);

    // OVERRIDDEN GETTERS
    const BaseONVBasis* get_fock_space() const override { return &fock_space; }
};


}  // namespace GQCP
