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


// #include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


#include <type_traits>


namespace GQCP {



// Forward declaration of EigenproblemEnvironment and EigenproblemSolver.
// This is needed because this step also requires the inclusion of EigenproblemSolver and EigenproblemEnvironment.

// class EigenproblemEnvironment {};
// class EigenproblemSolver {
// public:
//     static Algorithm<EigenproblemEnvironment> EigenproblemSolver::Dense(const size_t);
// };
// Algorithm<EigenproblemEnvironment> EigenproblemSolver::Dense(const size_t)




/**
 *  An iteration step that diagonalizes the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
 */
class SubspaceMatrixDiagonalization :
    public Step<EigenproblemEnvironment> {

private:
    size_t number_of_requested_eigenpairs;


public:

    /*
     *  CONSTRUCTORS
     */

    SubspaceMatrixDiagonalization(const size_t number_of_requested_eigenpairs = 1) :
        number_of_requested_eigenpairs (number_of_requested_eigenpairs)
    {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Calculate the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        // Diagonalize the subspace matrix and find the r (this->number_of_requested_eigenpairs) lowest eigenpairs
        // Lambda contains the requested number of eigenvalues, Z contains the corresponding eigenvectors
        // Z is a (subspace_dimension x number_of_requested_eigenpairs)- matrix

        // Use our dense diagonalization algorithm to find the number of requested eigenpairs
        const auto& S = environment.S;
        auto dense_environment = EigenproblemEnvironment::Dense(S);
        auto dense_diagonalizer = EigenproblemSolver::Dense(S.dimension());  // request all eigenpairs
        dense_diagonalizer.perform(dense_environment);


        environment.Lambda = dense_environment.eigenvalues.head(this->number_of_requested_eigenpairs);  // the (requested number of) eigenvalues of the subspace matrix S
        environment.eigenvalues = environment.Lambda;

        environment.Z = dense_environment.eigenvectors.topLeftCorner(S.cols(), this->number_of_requested_eigenpairs);  // the (requested number of) eigenvectors of the subspace matrix S
    }
};


}  // namespace GQCP
