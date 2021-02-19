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
#include "Mathematical/Optimization/Accelerator/DIIS.hpp"
#include "QCMethod/CC/CCDAmplitudesUpdate.hpp"
#include "QCMethod/CC/CCDIntermediatesUpdate.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"

#include <algorithm>


namespace GQCP {


/**
 *  An iteration step that accelerates the T2 amplitudes based on a DIIS accelerator.
 * 
 *  @tparam _Scalar              The scalar type used to represent the T2 amplitudes.
 */
template <typename _Scalar>
class T2DIIS:
    public Step<CCSDEnvironment<_Scalar>> {

public:
    // The scalar type used to represent the T2 amplitudes.
    using Scalar = _Scalar;

    // The type of environment that this iteration step can access.
    using Environment = CCSDEnvironment<Scalar>;


private:
    // The minimum number of T2 amplitudes that have to be in the subspace before enabling DIIS.
    size_t minimum_subspace_dimension;

    // The maximum number of T2 amplitues that can be handled by DIIS.
    size_t maximum_subspace_dimension;

    // The DIIS accelerator.
    DIIS<Scalar> diis;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param minimum_subspace_dimension       The minimum number of T2 amplitudes that have to be in the subspace before enabling DIIS.
     *  @param maximum_subspace_dimension       The maximum number of T2 amplitues that can be handled by DIIS.
     */
    T2DIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6) :
        minimum_subspace_dimension {minimum_subspace_dimension},
        maximum_subspace_dimension {maximum_subspace_dimension} {}


    /*
     *  MARK: Conforming to `Step`.
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Calculate the accelerated T2 amplitudes and place them in the environment by overwriting the previous T2 amplitudes.";
    }


    /**
     *  Calculate the accelerated T2 amplitudes and place them in the environment by overwriting the previous T2 amplitudes.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {

        // Don't do anything if the minimum number of T2 amplitude iterations isn't satisfied.
        if (environment.t2_amplitude_errors.size() < this->minimum_subspace_dimension) {
            return;
        }

        // Convert the deques in the environment to vectors that can be accepted by the DIIS accelerator. The total number of elements we can use in DIIS is either the maximum subspace dimension or the number of available error vectors.
        // TODO: Include the possibility for an x-iteration 'relaxation', i.e. not doing DIIS for x iterations long.
        const auto n = std::min(this->maximum_subspace_dimension, environment.t2_amplitude_errors.size());
        const std::vector<VectorX<Scalar>> error_vectors {environment.t2_amplitude_errors.end() - n, environment.t2_amplitude_errors.end()};  // The n-th last error vectors.
        const std::vector<T2Amplitudes<Scalar>> t2_amplitudes {environment.t2_amplitudes.end() - n, environment.t2_amplitudes.end()};         // The n-th last T2 amplitudes.


        // Calculate the accelerated T2 amplitudes and place them in the environment by overwriting the previous T2 amplitudes.
        const auto t2_amplitudes_accelerated = this->diis.accelerate(t2_amplitudes, error_vectors);

        environment.t2_amplitudes.pop_back();
        environment.t2_amplitudes.push_back(t2_amplitudes_accelerated);
    }
};


}  // namespace GQCP
