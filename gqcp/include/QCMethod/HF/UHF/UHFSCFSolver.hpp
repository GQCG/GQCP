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


#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/ConsecutiveIteratesNormConvergence.hpp"
#include "QCMethod/HF/RHF/UHFSCFEnvironment.hpp"


namespace GQCP {


/**
 *  A factory class that can construct UHF SCF solvers in an easy way.
 * 
 *  @tparam _Scalar             the scalar type that is used for the coefficient matrix/expansion coefficients
 */
template <typename _Scalar>
class UHFSCFSolver {
public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC STATIC METHODS
     */

    /**
     *  @param threshold                            the threshold that is used in comparing both the alpha and beta density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a plain UHF SCF solver that uses the combination of norm of the difference of two consecutive alpha and beta density matrices as a convergence criterion
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain RHF SCF solver
        StepCollection<RHFSCFEnvironment<Scalar>> plain_uhf_scf_cycle {};
        // plain_uhf_scf_cycle
        // .add(RHFDensityMatrixCalculation<Scalar>())
        // .add(RHFFockMatrixCalculation<Scalar>())
        // .add(RHFFockMatrixDiagonalization<Scalar>())
        // .add(RHFElectronicEnergyCalculation<Scalar>());

        // Create a convergence criterion on the norm of subsequent density matrices
        const auto density_matrix_extractor = [](const RHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices; };

        using ConvergenceType = ConsecutiveIteratesNormConvergence<OneRDM<Scalar>, RHFSCFEnvironment<Scalar>>;
        const ConvergenceType convergence_criterion {threshold, density_matrix_extractor};

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(plain_rhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
