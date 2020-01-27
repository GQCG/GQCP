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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Accelerator/ConstantDamper.hpp"
#include "QCMethod/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/RHF/RHF.hpp"


namespace GQCP {


/**
 *  An iteration step that accelerates the density matrix (expressed in the scalar/AO basis) based on a constant damping accelerator.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class RHFDensityMatrixDamper :
    public Step<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


private:
    ConstantDamper damper;  // the damping accelerator


public:

    /*
     *  CONSTRUCTORS 
     */

    /**
     *  @param alpha            the damping factor
     */
    RHFDensityMatrixDamper(const double alpha) : 
        damper (alpha)
    {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Replace the most recent density matrix with an accelerated one.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        if (environment.density_matrices.size() < 2) {
            return;  // no acceleration is possible
        }

        // Get the two most recent density matrices and produce an accelerated density matrix
        const auto second_to_last_it = environment.density_matrices.end() - 2;
        const auto& D_previous = *second_to_last_it;  // dereference the iterator
        const auto D_current = environment.density_matrices.back();

        const auto D_accelerated = this->damper.accelerate(D_current, D_previous);
        environment.density_matrices.pop_back();  // we will replace the most recent density matrix with the accelerated one, so remove the most recent one
        environment.density_matrices.push_back(D_accelerated);
    }
};


}  // namespace GQCP
