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


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/CC/SinglesAmplitudes.hpp"


namespace GQCP {


/**
 *  The coupled-cluster Lambda-1 amplitudes l_i^a.
 * 
 *  @param _Scalar:             The scalar type of one of the amplitudes.
 */
template <typename _Scalar>
class L1Amplitudes:
    public SinglesAmplitudes<_Scalar, L1Amplitudes<_Scalar>> {

public:
    // The scalar type of one of the amplitudes.
    using Scalar = _Scalar;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SingleAmplitudes`' constructors.
    using SingleAmplitudes<Scalar, L1Amplitudes<Scalar>>::SingleAmplitudes;
};


}  // namespace GQCP
