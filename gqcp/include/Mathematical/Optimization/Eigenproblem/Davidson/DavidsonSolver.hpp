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
#include "Mathematical/Algorithm/StepCollection.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/CorrectionVectorCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/GuessVectorUpdate.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/MatrixVectorProductCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/ResidualVectorCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/ResidualVectorConvergence.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceMatrixCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceMatrixDiagonalization.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceUpdate.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


namespace GQCP {
namespace EigenproblemSolver {


/**
 *  @param number_of_requested_eigenpairs       the number of solutions the Davidson solver should find
 *  @param maximum_subspace_dimension           the maximum dimension of the subspace before collapsing
 *  @param convergence_threshold                the threshold that is used in determining the norm on the residuals, which determines convergence
 *  @param correction_threshold                 the threshold used in solving the (approximated) residue correction equation
 *  @param maximum_number_of_iterations         the maximum number of iterations the algorithm may perform
 *  @param inclusion_threshold                  the threshold on the norm used for determining if a new projected correction vector should be added to the subspace
 * 
 *  @return an iterative algorithm that can find the lowest n eigenvectors of a matrix using Davidson's algorithm
 */
IterativeAlgorithm<EigenproblemEnvironment> Davidson(const size_t number_of_requested_eigenpairs = 1, const size_t maximum_subspace_dimension = 15, const double convergence_threshold = 1.0e-08, double correction_threshold = 1.0e-12, const size_t maximum_number_of_iterations = 128, const double inclusion_threshold = 1.0e-03) {

    // Create the iteration cycle that effectively 'defines' our Davidson solver
    StepCollection<EigenproblemEnvironment> davidson_cycle {};

    davidson_cycle
        .add(MatrixVectorProductCalculation())
        .add(SubspaceMatrixCalculation())
        .add(SubspaceMatrixDiagonalization(number_of_requested_eigenpairs))
        .add(GuessVectorUpdate())
        .add(ResidualVectorCalculation(number_of_requested_eigenpairs))
        .add(CorrectionVectorCalculation(number_of_requested_eigenpairs, correction_threshold))  // this solves the residual equations
        .add(SubspaceUpdate(maximum_subspace_dimension, inclusion_threshold));

    // Create a convergence criterion on the norm of the residual vectors
    const ResidualVectorConvergence<EigenproblemEnvironment> convergence_criterion {convergence_threshold};

    return IterativeAlgorithm<EigenproblemEnvironment>(davidson_cycle, convergence_criterion, maximum_number_of_iterations);
}


}  // namespace EigenproblemSolver
}  // namespace GQCP
