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


#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  An iteration step that diagonalizes the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V.
 */
class SubspaceMatrixDiagonalization:
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
    SubspaceMatrixDiagonalization(const size_t number_of_requested_eigenpairs = 1) :
        number_of_requested_eigenpairs {number_of_requested_eigenpairs} {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Diagonalize the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V, and write its eigenvalues and eigenvectors to the environment.";
    }


    /**
     *  Diagonalize the subspace matrix, i.e. the projection of the matrix A onto the subspace spanned by the vectors in V, and write its eigenvalues and eigenvectors to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(EigenproblemEnvironment& environment) override {

        // Diagonalize the subspace matrix and find the r (this->number_of_requested_eigenpairs) lowest eigenpairs
        // Lambda contains the requested number of eigenvalues, Z contains the corresponding eigenvectors
        // Z is a (subspace_dimension x number_of_requested_eigenpairs)- matrix

        // Use our own dense diagonalization algorithm to find the number of requested eigenpairs
        const auto& S = environment.S;
        auto dense_environment = EigenproblemEnvironment::Dense(S);
        auto dense_diagonalizer = EigenproblemSolver::Dense();
        dense_diagonalizer.perform(dense_environment);


        environment.Lambda = dense_environment.eigenvalues.head(this->number_of_requested_eigenpairs);  // the (requested number of) eigenvalues of the subspace matrix S
        environment.eigenvalues = environment.Lambda;

        environment.Z = dense_environment.eigenvectors.topLeftCorner(S.cols(), this->number_of_requested_eigenpairs);  // the (requested number of) eigenvectors of the subspace matrix S
    }
};


}  // namespace GQCP
