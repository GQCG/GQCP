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
#pragma once


#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreCI.hpp"


namespace GQCP {


/**
 *  A class capable of generating the matrix representation of the frozen core FCI Hamiltonian
 */
class FrozenCoreFCI : public FrozenCoreCI {
private:
    SpinResolvedFrozenONVBasis onv_basis;  // the frozen spin-resolved ONV

public:
    // CONSTRUCTORS
    /**
     *  @param onv_basis       the frozen spin-resolved ONV
     */
    explicit FrozenCoreFCI(const SpinResolvedFrozenONVBasis& onv_basis);


    // OVERRIDDEN GETTERS
    const BaseONVBasis* get_fock_space() const override { return &onv_basis; }
};


}  // namespace GQCP
