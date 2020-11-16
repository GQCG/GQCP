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
#include "QCMethod/HF/GHF/GHFSCFEnvironment.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  An iteration step that solves the generalized eigenvalue problem for the current scalar/AO basis Fock matrix for the coefficient matrix.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class GHFFockMatrixDiagonalization:
    public Step<GHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = GHFSCFEnvironment<Scalar>;


public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Solve the generalized eigenvalue problem for the most recent scalar/AO Fock matrix. Add the associated coefficient matrix and orbital energies to the environment.";
    }


    /**
     *  Solve the generalized eigenvalue problem for the most recent scalar/AO Fock matrix. Add the associated coefficient matrix and orbital energies to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        const auto& F = environment.fock_matrices.back().parameters();  // the most recent scalar/AO basis Fock matrix

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver {F, environment.S.parameters()};
        const GTransformationMatrix<Scalar>& C = generalized_eigensolver.eigenvectors();
        const auto& orbital_energies = generalized_eigensolver.eigenvalues();

        environment.coefficient_matrices.push_back(C);
        environment.orbital_energies.push_back(orbital_energies);
    }
};


}  // namespace GQCP
