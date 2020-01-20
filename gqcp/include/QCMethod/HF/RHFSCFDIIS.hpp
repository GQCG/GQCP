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


#include "Mathematical/Optimization/Accelerators/DIIS.hpp"


namespace GQCP {


/**
 *  A DIIS accelerator for RHF SCF.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class RHFSCFDIIS :
    public DIIS<ScalarSQOneElectronOperator<_ExpansionScalar>, ScalarSQOneElectronOperator<_ExpansionScalar>, ExpansionScalar> {

public:
    using ExpansionScalar = _ExpansionScalar;



public:

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @param e1                   the first error
     *  @param e2                   the second error
     * 
     *  @return the scalar product between the two given errors
     */
    virtual ExpansionScalar calculateErrorScalarProduct(const Error& e1, const Error& e2);
};



}  // namespace GQCP