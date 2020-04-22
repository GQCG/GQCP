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


/**
 *  A step that calculates the residual vectors from the new guesses for the eigenvectors.
 */
class ResidualVectorCalculation:
    public Step<EigenproblemEnvironment> {

private:
    size_t number_of_requested_eigenpairs;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param number_of_requested_eigenpairs       the number of solutions the Davidson solver should find
     */
    ResidualVectorCalculation(const size_t number_of_requested_eigenpairs = 1) :
        number_of_requested_eigenpairs{number_of_requested_eigenpairs} {}


    /**
     *  Calculate the residual vectors from the new guesses for the eigenvectors.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        const auto dim = environment.dimension;

        const auto& VA = environment.VA;          // VA = A * V (implicitly calculated through the matrix-vector product)
        const auto& Z = environment.Z;            // the (requested number of) eigenvectors of the subspace matrix S
        const auto& Lambda = environment.Lambda;  // the (requested number of) eigenvalues of the subspace matrix S
        const auto& X = environment.X;            // contains the new guesses for the eigenvectors (as a linear combination of the current subspace V)

        // Calculate the residual vectors: r_i = VA * z_i - Lambda * x_i
        environment.R = MatrixX<double>::Zero(dim, this->number_of_requested_eigenpairs);
        for (size_t column_index = 0; column_index < this->number_of_requested_eigenpairs; column_index++) {
            environment.R.col(column_index) = VA * Z.col(column_index) - Lambda(column_index) * X.col(column_index);
        }
    }
};


}  // namespace GQCP
