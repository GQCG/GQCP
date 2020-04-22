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


#include "ONVBasis/SpinUnresolvedONV.hpp"


namespace GQCP {


/**
 *  A spin-resolved ONV.
 */
class SpinResolvedONV {
private:
    SpinUnresolvedONV onv_alpha;  // the ONV that descripes the occupations of the alpha spinors
    SpinUnresolvedONV onv_beta;   // the ONV that descripes the occupations of the beta spinors

public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param onv_alpha                the ONV that descripes the occupations of the alpha spinors
     *  @param onv_beta                 the ONV that descripes the occupations of the beta spinors
     */
    SpinResolvedONV(const SpinUnresolvedONV& onv_alpha, const SpinUnresolvedONV& onv_beta) :
        onv_alpha {onv_alpha},
        onv_beta {onv_beta} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the ONV that descripes the occupations of the alpha spinors
     */
    const SpinUnresolvedONV& alphaONV() const { return this->onv_alpha; }

    /**
     *  @return the ONV that descripes the occupations of the beta spinors
     */
    const SpinUnresolvedONV& betaONV() const { return this->onv_beta; }
};


}  // namespace GQCP
