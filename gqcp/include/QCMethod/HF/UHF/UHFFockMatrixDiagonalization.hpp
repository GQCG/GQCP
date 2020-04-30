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
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  An iteration step that solves the generalized eigenvalue problem for the current scalar/AO basis Fock matrix for the coefficient matrix.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class UHFFockMatrixDiagonalization:
    public Step<UHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = UHFSCFEnvironment<Scalar>;


public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  Solve the generalized eigenvalue problem for the most recent scalar/AO Fock matrices. Add the associated coefficient matrices and orbital energies to the environment.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        const auto& F_alpha = environment.fock_matrices_alpha.back();  // the most recent scalar/AO basis alpha Fock matrix
        const auto& F_beta = environment.fock_matrices_beta.back();    // the most recent scalar/AO basis beta Fock matrices

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver_alpha {F_alpha, environment.S};
        const TransformationMatrix<Scalar>& C_alpha = generalized_eigensolver_alpha.eigenvectors();
        environment.coefficient_matrices_alpha.push_back(C_alpha);

        const auto& orbital_energies_alpha = generalized_eigensolver_alpha.eigenvalues();
        environment.orbital_energies_alpha.push_back(orbital_energies_alpha);


        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver_beta {F_beta, environment.S};
        const TransformationMatrix<Scalar>& C_beta = generalized_eigensolver_beta.eigenvectors();
        environment.coefficient_matrices_beta.push_back(C_beta);

        const auto& orbital_energies_beta = generalized_eigensolver_beta.eigenvalues();
        environment.orbital_energies_beta.push_back(orbital_energies_beta);
    }
};


}  // namespace GQCP
