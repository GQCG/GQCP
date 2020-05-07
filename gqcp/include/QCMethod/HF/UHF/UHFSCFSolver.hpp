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


#include "Mathematical/Algorithm/CompoundConvergenceCriterion.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/ConsecutiveIteratesNormConvergence.hpp"
#include "QCMethod/HF/UHF/UHFDensityMatrixCalculation.hpp"
#include "QCMethod/HF/UHF/UHFElectronicEnergyCalculation.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixCalculation.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"


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
    static IterativeAlgorithm<UHFSCFEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain UHF SCF solver
        StepCollection<UHFSCFEnvironment<Scalar>> plain_uhf_scf_cycle {};
        plain_uhf_scf_cycle
            .add(UHFDensityMatrixCalculation<Scalar>())
            .add(UHFFockMatrixCalculation<Scalar>())
            .add(UHFFockMatrixDiagonalization<Scalar>())
            .add(UHFElectronicEnergyCalculation<Scalar>());

        // Create a compound convergence criterion on the norm of subsequent alpha- and beta-density matrices
        using SingleConvergenceType = ConsecutiveIteratesNormConvergence<OneRDM<Scalar>, UHFSCFEnvironment<Scalar>>;

        const auto density_matrix_alpha_extractor = [](const UHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices_alpha; };
        const SingleConvergenceType convergence_criterion_alpha {threshold, density_matrix_alpha_extractor, "the UHF alpha-density matrix in AO basis"};

        const auto density_matrix_beta_extractor = [](const UHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices_beta; };
        const SingleConvergenceType convergence_criterion_beta {threshold, density_matrix_beta_extractor, "the UHF beta-density matrix in AO basis"};

        const CompoundConvergenceCriterion<UHFSCFEnvironment<Scalar>> convergence_criterion {convergence_criterion_alpha, convergence_criterion_beta};


        // Put together the pieces of the algorithm.
        return IterativeAlgorithm<UHFSCFEnvironment<Scalar>>(plain_uhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
