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
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreCI.hpp"


namespace GQCP {


/**
 *  A class capable of generating the matrix representation of the frozen core FCI Hamiltonian
 */
class FrozenCoreFCI: public FrozenCoreCI {
private:
    SpinResolvedFrozenONVBasis onv_basis;  // the frozen spin-resolved ONV


public:
    // CONSTRUCTORS

    /**
     *  @param onv_basis       the frozen spin-resolved ONV
     */
    explicit FrozenCoreFCI(const SpinResolvedFrozenONVBasis& onv_basis);


    // PUBLIC METHODS

    /**
     *  @return the ONV basis that is associated to this HamiltonianBuilder
     */
    const BaseONVBasis* onvBasis() const override { return &onv_basis; }
};


}  // namespace GQCP
