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


#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


namespace GQCP {


class GeminalCoefficientsInterface {
public:
    // DESTRUCTOR
    virtual ~GeminalCoefficientsInterface() = default;


    // PUBLIC METHODS
    /**
     *  @param onv      the doubly-occupied (spin-resolved) ONV that is being projected on
     *
     *  @return the overlap of the APIG-like wave function with the given ONV, i.e. the projection of the APIG wave function onto that ONV
     */
    virtual double overlap(const SpinUnresolvedONV& onv) const = 0;

    /**
     *  @param onv_basis       the seniority-zero spin-resolved ONV basis the wave function should live in
     *
     *  @return the wave function expansion corresponding to the geminal coefficients
     */
    LinearExpansion<SeniorityZeroONVBasis> toLinearExpansion(const SeniorityZeroONVBasis& onv_basis) const;
};


}  // namespace GQCP
