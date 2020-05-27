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


#include "Mathematical/Algorithm/Step.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"


namespace GQCP {


/**
 *  An iteration step that calculates the current CCD intermediates as described in Stanton1991.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the amplitudes
 */
template <typename _Scalar>
class CCDIntermediatesUpdate:
    public Step<CCSDEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = CCSDEnvironment<Scalar>;


public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the current CCD intermediates as described in Stanton1991.";
    }

    /**
     *  Calculate the current CCD intermediates as described in Stanton1991.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        // Extract and prepare some variables.
        const auto& f = environment.f;
        const auto& V_A = environment.V_A;
        const auto& t2 = environment.t2_amplitudes.back();

        // Calculate the other CCD intermediates and push them to the environment.
        environment.F1 = QCModel::CCD<Scalar>::calculateF1(f, V_A, t2);
        environment.F2 = QCModel::CCD<Scalar>::calculateF2(f, V_A, t2);

        environment.W1 = QCModel::CCD<Scalar>::calculateW1(V_A, t2);
        environment.W2 = QCModel::CCD<Scalar>::calculateW2(V_A, t2);
        environment.W3 = QCModel::CCD<Scalar>::calculateW3(V_A, t2);
    }
};


}  // namespace GQCP
