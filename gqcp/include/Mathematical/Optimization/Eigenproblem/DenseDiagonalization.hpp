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

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A step that performs a dense diagonalization.
 */
class DenseDiagonalization:
    public Step<EigenproblemEnvironment> {

private:
    size_t number_of_requested_eigenpairs;


public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Read the matrix from the environment, diagonalize it and write the number of requested eigenpairs to it
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(Environment& environment) override {

        const auto& A = environment.A;                                  // a self-adjoint matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(A);  // this solves the eigenvalue problem

        // Write the eigenvalues and eigenvectors to the environment
        environment.eigenvalues = eigensolver.eigenvalues();
        environment.eigenvectors = eigensolver.eigenvectors();
    }
};


}  // namespace GQCP
