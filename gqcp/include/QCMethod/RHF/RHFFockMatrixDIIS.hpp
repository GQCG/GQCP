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


#include "Mathematical/Algorithm/IterationStep.hpp"
#include "Mathematical/Optimization/Accelerator/DIIS.hpp"
#include "QCMethod/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/RHF/RHF.hpp"
#include "QCMethod/RHF/RHFFockMatrixDiagonalization.hpp"
#include "QCMethod/RHF/RHFFockMatrixDiagonalization.hpp"

#include <algorithm>


namespace GQCP {


/**
 *  An iteration step that accelerates the density matrix (expressed in the scalar/AO basis) based on a DIIS accelerator.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class RHFFockMatrixDIIS :
    public IterationStep<RHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = RHFSCFEnvironment<Scalar>;


private:
    size_t minimum_subspace_dimension;  // the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
    size_t maximum_subspace_dimension;  // the maximum number of Fock matrices that can be handled by DIIS

    DIIS<Scalar> diis;  // the DIIS accelerator


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param minimum_subspace_dimension       the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
     *  @param maximum_subspace_dimension       the maximum number of Fock matrices that can be handled by DIIS
     */
    RHFFockMatrixDIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6) :
        minimum_subspace_dimension (minimum_subspace_dimension),
        maximum_subspace_dimension (maximum_subspace_dimension)
    {}


    /*
     *  PUBLIC METHODS
     */

    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the accelerated Fock matrix, and perform a diagonalization step on it.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        if (environment.error_vectors.size() < this->minimum_subspace_dimension) {

            // No acceleration is possible, so calculate the regular Fock matrix and diagonalize it.
            RHFFockMatrixCalculation<Scalar>().execute(environment);
            RHFFockMatrixDiagonalization<Scalar>().execute(environment);
            return;
        }

        // Convert the deques in the environment to vectors that can be accepted by the DIIS accelerator. The total number of elements we can use in DIIS is either the maximum subspace dimension or the number of available error matrices.
        const auto n = std::min(this->maximum_subspace_dimension, environment.error_vectors.size());
        const std::vector<VectorX<Scalar>> error_vectors (environment.error_vectors.end() - n, environment.error_vectors.end());  // the n-th last error vectors
        const std::vector<QCMatrix<Scalar>> fock_matrices (environment.fock_matrices.end() - n, environment.fock_matrices.end());  // the n-th last Fock matrices

        // Calculate the accelerated Fock matrix and do a diagonalization step on it
        const auto F_accelerated = this->diis.accelerate(fock_matrices, error_vectors);

        environment.fock_matrices.push_back(F_accelerated);  // the diagonalization step can only read from the environment
        RHFFockMatrixDiagonalization<Scalar>().execute(environment);
        environment.fock_matrices.pop_back();  // the accelerated/extrapolated Fock matrix should not be used in further extrapolation steps, as it is not created from a density matrix
    }
};


}  // namespace GQCP
