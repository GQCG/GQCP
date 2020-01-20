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


#include "Mathematical/Algorithm/ConvergenceCriterion.hpp"
#include "QCMethod/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/RHF/RHF.hpp"


namespace GQCP {


/**
 *  A convergence criterion on the norm of subsequent RHF density matrices.
 * 
 *  @param _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class RHFDensityMatrixConvergenceCriterion :
    public ConvergenceCriterion<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


private:
    double threshold;  // the threshold that is used in comparing the density matrices


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param threshold                    the threshold that is used in comparing the density matrices
     */
    RHFDensityMatrixConvergenceCriterion(const double threshold = 1.0e-08) :
        threshold (threshold)
    {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @param environment              the environment that acts as a sort of calculation space
     * 
     *  @return if the difference of the two most recent density matrices has a zero norm, within the tolerance
     */
    bool isFulfilled(Environment& environment) override {

        if (environment.density_matrices.size() < 2) {
            return false;
        }

        // Get the two most recent density matrices and compare the norm of their difference
        const auto second_to_last_it = environment.density_matrices.end() - 1;
        const auto& D_previous = *second_to_last_it;  // dereference the iterator
        const auto D_current = environment.density_matrices.back();

        return (std::abs((D_current - D_previous).norm()) < this->threshold);
    }
};


}  // namespace GQCP
