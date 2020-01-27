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

#include <Eigen/Dense>


namespace GQCP {


/**
 * 
 */
class DenseDiagonalization :
    public Step<EigenproblemEnvironment> {

private:
    size_t number_of_requested_eigenpairs;


public:

    /*
     *  CONSTRUCTORS
     */

    DenseDiagonalization(const size_t number_of_requested_eigenpairs = 1) :
        number_of_requested_eigenpairs (number_of_requested_eigenpairs)
    {}


    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Read the matrix from the environment, diagonalize it and write the number of requested eigenpairs to it
     * 
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(Environment& environment) override {

        const auto& A = environment.A;  // a self-adjoint matrix
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver (A);  // this solves the eigenvalue problem

        // Write the eigenpairs back to the environment
        environment.eigenpairs.reserve(this->number_of_requested_eigenpairs);
        for (size_t i = 0; i < this->number_of_requested_eigenpairs; i++) {
            const auto& eigenvalue = eigensolver.eigenvalues()(i);
            const auto& eigenvector = eigensolver.eigenvectors().col(i);

            environment.eigenpairs.emplace_back(eigenvalue, eigenvector);
        }
    }
};


}  // namespace GQCP
