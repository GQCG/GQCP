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


#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "QCMethod/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/RHF/RHFDensityMatrixCalculation.hpp"
#include "QCMethod/RHF/RHFDensityMatrixConvergenceCriterion.hpp"
#include "QCMethod/RHF/RHFDensityMatrixDamper.hpp"
#include "QCMethod/RHF/RHFErrorCalculation.hpp"
#include "QCMethod/RHF/RHFElectronicEnergyCalculation.hpp"
#include "QCMethod/RHF/RHFFockMatrixCalculation.hpp"
#include "QCMethod/RHF/RHFFockMatrixDiagonalization.hpp"
#include "QCMethod/RHF/RHFFockMatrixDIIS.hpp"


namespace GQCP {


/**
 *  A factory class that can construct RHF SCF solvers in an easy way.
 * 
 *  @tparam _Scalar             the scalar type that is used for the coefficient matrix/expansion coefficients
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
     *  @param threshold                            the threshold that is used in comparing the density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a plain RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> Plain(const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a plain RHF SCF solver
        IterationCycle<RHFSCFEnvironment<Scalar>> plain_rhf_scf_cycle {};
        plain_rhf_scf_cycle.add(RHFDensityMatrixCalculation<Scalar>())
                           .add(RHFFockMatrixCalculation<Scalar>())
                           .add(RHFFockMatrixDiagonalization<Scalar>())
                           .add(RHFElectronicEnergyCalculation<Scalar>());

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(plain_rhf_scf_cycle, RHFDensityMatrixConvergenceCriterion<Scalar>(threshold), maximum_number_of_iterations);
    }


    /**
     *  @param alpha                                the damping factor
     *  @param threshold                            the threshold that is used in comparing the density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a density-damped RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> DensityDamped(const double alpha, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a damped RHF SCF solver
        IterationCycle<RHFSCFEnvironment<Scalar>> damped_rhf_scf_cycle {};
        damped_rhf_scf_cycle.add(RHFDensityMatrixCalculation<Scalar>())
                            .add(RHFDensityMatrixDamper<Scalar>(alpha))
                            .add(RHFFockMatrixCalculation<Scalar>())
                            .add(RHFFockMatrixDiagonalization<Scalar>())
                            .add(RHFElectronicEnergyCalculation<Scalar>());

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(damped_rhf_scf_cycle, RHFDensityMatrixConvergenceCriterion<Scalar>(threshold), maximum_number_of_iterations);
    }


    /**
     *  @param threshold                            the threshold that is used in comparing the density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
     * 
     *  @return a DIIS RHF SCF solver that uses the norm of the difference of two consecutive density matrices as a convergence criterion.
     */
    static IterativeAlgorithm<RHFSCFEnvironment<Scalar>> DIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {

        // Create the iteration cycle that effectively 'defines' a DIIS RHF SCF solver
        IterationCycle<RHFSCFEnvironment<Scalar>> diis_rhf_scf_cycle {};
        diis_rhf_scf_cycle.add(RHFDensityMatrixCalculation<Scalar>())
                          .add(RHFFockMatrixCalculation<Scalar>())
                          .add(RHFErrorCalculation<Scalar>())
                          .add(RHFFockMatrixDIIS<Scalar>(minimum_subspace_dimension, maximum_subspace_dimension))  // this also calculates the next coefficient matrix
                          .add(RHFElectronicEnergyCalculation<Scalar>());

        return IterativeAlgorithm<RHFSCFEnvironment<Scalar>>(diis_rhf_scf_cycle, RHFDensityMatrixConvergenceCriterion<Scalar>(threshold), maximum_number_of_iterations);
    } 
};


}  // namespace GQCP
