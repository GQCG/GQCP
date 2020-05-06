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


#include "ONVBasis/SpinUnresolvedFrozenONVBasis.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreCI.hpp"


namespace GQCP {


/**
 *  A class capable of generating the matrix representation of the frozen core DOCI Hamiltonian
 */
class FrozenCoreDOCI: public FrozenCoreCI {
private:
    SpinUnresolvedFrozenONVBasis onv_basis;  // both the alpha and beta frozen ONV basis

public:
    // CONSTRUCTORS
    /**
     *  @param onv_basis       the frozen ONV basis (identical for alpha and beta)
     */
    explicit FrozenCoreDOCI(const SpinUnresolvedFrozenONVBasis& onv_basis);


    // OVERRIDDEN GETTERS
    const BaseONVBasis* get_fock_space() const override { return &onv_basis; }
};


}  // namespace GQCP
