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
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


namespace GQCP {


class CorrectionVectorCalculation : 
    public Step<EigenproblemEnvironment> {

private:
    size_t number_of_requested_eigenpairs;
    double correction_threshold;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param number_of_requestested_eigenpairs            the number of eigenpairs that should be found by the algorithm
     *  @param correction_threshold                         the threshold used in solving the (approximated) residue correction equation
     */
    CorrectionVectorCalculation(const size_t number_of_requested_eigenpairs = 1, const double correction_threshold = 1.0e-12) :
        number_of_requested_eigenpairs (number_of_requested_eigenpairs),
        correction_threshold (correction_threshold)
    {}


    /**
     *  Calculate the correction vectors by solving the residual equations.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        const auto dim = environment.dimension;

        const auto& diagonal = environment.diagonal;  // the diagonal of the matrix
        const auto& Lambda = environment.Lambda;  // the (requested number of) eigenvalues of the subspace matrix S

        const auto& VA = environment.VA;  // VA = A * V (implicitly calculated through the matrix-vector product)
        const auto& Z = environment.Z;  // the (requested number of) eigenvectors of the subspace matrix S
        const auto& X = environment.X;  // contains the new guesses for the eigenvectors (as a linear combination of the current subspace V)
        const auto& R = environment.R;  // the residual vectors

        // Solve the residual equations to find the correction vectors.
        // The implementation of these equations is adapted from Klaas Gunst's DOCI code (https://github.com/klgunst/doci)
        environment.Delta = MatrixX<double>::Zero(dim, this->number_of_requested_eigenpairs);
        for (size_t column_index = 0; column_index < this->number_of_requested_eigenpairs; column_index++) {

            VectorX<double> denominator = diagonal - VectorX<double>::Constant(dim, Lambda(column_index));

            // If the denominator is large enough, the correction vector is the residual vector dividided by the denominator.
            // If it isn't, the correction vector is the residual vector divided by the threshold.
            environment.Delta.col(column_index) = (denominator.array().abs() > this->correction_threshold).select( 
                R.col(column_index).array() / denominator.array().abs(), 
                R.col(column_index) / this->correction_threshold
            );
            environment.Delta.col(column_index).normalize();
        }
    }
};


}  // namespace GQCP
