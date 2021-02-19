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
 *  An iteration step that calculates the current T2 amplitude error.
 * 
 *  @tparam _Scalar              The scalar type used to represent the T2 amplitudes.
 */
template <typename _Scalar>
class T2ErrorCalculation:
    public Step<CCSDEnvironment<_Scalar>> {

public:
    // The scalar type used to represent the T2 amplitudes.
    using Scalar = _Scalar;

    // The environment related to this step.
    using Environment = CCSDEnvironment<Scalar>;


public:
    /*
     *  MARK: Conforming to `Step`
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Calculate the current T2 error vector and add it to the environment.";
    }


    /**
     *  Calculate the current T2 error vector and add it to the environment.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {

        // Read the last two T2 amplitudes iterations and calculate the error as their difference.
        const auto second_to_last_it = environment.t2_amplitudes.end() - 2;
        const auto& T2_previous = *second_to_last_it;  // Dereference the iterator.

        const auto& T2_current = environment.t2_amplitudes.back();

        // Calculate the current T2 error vector and add it to the environment (as a vector).
        const auto t2_error = T2_current - T2_previous;
        environment.t2_amplitude_errors.push_back(t2_error.asImplicitRankFourTensorSlice().asMatrix().pairWiseReduced());
    }
};


}  // namespace GQCP
