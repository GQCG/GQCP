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
#include "QCMethod/HF/RHF/RHFDensityMatrixCalculation.hpp"
#include "QCMethod/HF/RHF/RHFDensityMatrixDamper.hpp"
#include "QCMethod/HF/RHF/RHFElectronicEnergyCalculation.hpp"
#include "QCMethod/HF/RHF/RHFErrorCalculation.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixCalculation.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixDIIS.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"


namespace GQCP {


/**
 *  A factory class that can construct RHF SCF solvers in an easy way.
 * 
 *  @tparam _Scalar             The scalar type that is used for the coefficient matrix/expansion coefficients: real or complex.
 */
template <typename _Scalar>
class RHFSCFSolver {
public:
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC STATIC METHODS
     */

    /**
     *  @param alpha                                The damping factor.
     *  @param threshold                            The threshold that is used in comparing the density matrices.
     *  @param maximum_number_of_iterations         The maximum number of iterations the algorithm may perform.
     * 
     *  @return A density-damped RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> DensityDamped(const double alpha, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a damped RHF SCF solver.
        StepCollection<RHFSCFEnvironment<Scalar>> damped_rhf_scf_cycle {};
        damped_rhf_scf_cycle
            .add(RHFDensityMatrixCalculation<Scalar>())
            .add(RHFDensityMatrixDamper<Scalar>(alpha))
            .add(RHFFockMatrixCalculation<Scalar>())
            .add(RHFFockMatrixDiagonalization<Scalar>())
            .add(RHFElectronicEnergyCalculation<Scalar>());

        // Create a convergence criterion on the norm of subsequent density matrices.
        const std::function<std::deque<Orbital1DM<Scalar>>(const RHFSCFEnvironment<Scalar>&)> density_matrix_extractor = [](const RHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices; };

        using ConvergenceType = ConsecutiveIteratesNormConvergence<Orbital1DM<Scalar>, RHFSCFEnvironment<Scalar>>;
        const ConvergenceType convergence_criterion {threshold, density_matrix_extractor, "the RHF density matrix in AO basis"};

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(damped_rhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
    }


    /**
     *  @param minimum_subspace_dimension           The minimum number of Fock matrices that have to be in the subspace before enabling DIIS.
     *  @param maximum_subspace_dimension           The maximum number of Fock matrices that can be handled by DIIS.
     *  @param threshold                            The threshold that is used in comparing the density matrices.
     *  @param maximum_number_of_iterations         The maximum number of iterations the algorithm may perform.
     * 
     *  @return A DIIS RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> DIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a DIIS RHF SCF solver.
        StepCollection<RHFSCFEnvironment<Scalar>> diis_rhf_scf_cycle {};
        diis_rhf_scf_cycle
            .add(RHFDensityMatrixCalculation<Scalar>())
            .add(RHFFockMatrixCalculation<Scalar>())
            .add(RHFErrorCalculation<Scalar>())
            .add(RHFFockMatrixDIIS<Scalar>(minimum_subspace_dimension, maximum_subspace_dimension))  // This also calculates the next coefficient matrix.
            .add(RHFElectronicEnergyCalculation<Scalar>());

        // Create a convergence criterion on the norm of subsequent density matrices.
        const std::function<std::deque<Orbital1DM<Scalar>>(const RHFSCFEnvironment<Scalar>&)> density_matrix_extractor = [](const RHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices; };

        using ConvergenceType = ConsecutiveIteratesNormConvergence<Orbital1DM<Scalar>, RHFSCFEnvironment<Scalar>>;
        const ConvergenceType convergence_criterion {threshold, density_matrix_extractor, "the RHF density matrix in AO basis"};

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(diis_rhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
    }


    /**
     *  @param threshold                            The threshold that is used in comparing the density matrices.
     *  @param maximum_number_of_iterations         The maximum number of iterations the algorithm may perform.
     * 
     *  @return A plain RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain RHF SCF solver.
        StepCollection<RHFSCFEnvironment<Scalar>> plain_rhf_scf_cycle {};
        plain_rhf_scf_cycle
            .add(RHFDensityMatrixCalculation<Scalar>())
            .add(RHFFockMatrixCalculation<Scalar>())
            .add(RHFFockMatrixDiagonalization<Scalar>())
            .add(RHFElectronicEnergyCalculation<Scalar>());

        // Create a convergence criterion on the norm of subsequent density matrices.
        const auto density_matrix_extractor = [](const RHFSCFEnvironment<Scalar>& environment) { return environment.density_matrices; };

        using ConvergenceType = ConsecutiveIteratesNormConvergence<Orbital1DM<Scalar>, RHFSCFEnvironment<Scalar>>;
        const ConvergenceType convergence_criterion {threshold, density_matrix_extractor, "the RHF density matrix in AO basis"};

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(plain_rhf_scf_cycle, convergence_criterion, maximum_number_of_iterations);
    }
};


}  // namespace GQCP
