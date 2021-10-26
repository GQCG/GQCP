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
#include "Mathematical/Optimization/Eigenproblem/GeneralizedEigenproblemEnvironment.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  A step that performs a dense diagonalization.
 *
 *  @tparam _Scalar         The scalar type of the matrix elements: real or complex.
 */
template <typename _Scalar>
class GeneralizedDenseDiagonalization:
    public Step<GeneralizedEigenproblemEnvironment<_Scalar>> {
public:
    // The scalar type of the matrix elements: real or complex.
    using Scalar = _Scalar;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Read the matrix from the environment, diagonalize it and write the number of requested eigenpairs to it.";
    }


    /**
     *  Read the matrix from the environment, diagonalize it and write the number of requested eigenpairs to it.
     *
     *  @param environment              the environment that this step can read from and write to
     */
    void execute(GeneralizedEigenproblemEnvironment<Scalar>& environment) override {

        // Solve the eigenvalue problem using Eigen's routines.
        const auto& A = environment.A;
        const auto& S = environment.S;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> eigensolver {A, S};

        // Write the eigenvalues and eigenvectors to the environment.
        environment.eigenvalues = eigensolver.eigenvalues();
        environment.eigenvectors = eigensolver.eigenvectors();
    }
};


}  // namespace GQCP
