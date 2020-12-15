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
 *  @tparam _Scalar              The scalar type used to represent the expansion coefficient/elements of the transformation matrix: real or complex.
 */
template <typename _Scalar>
class UHFFockMatrixDiagonalization:
    public Step<UHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = UHFSCFEnvironment<Scalar>;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return A textual description of this algorithmic step.
     */
    std::string description() const override {
        return "Solve the generalized eigenvalue problem for the most recent scalar/AO Fock matrices. Add the associated coefficient matrices and orbital energies to the environment.";
    }


    /**
     *  Solve the generalized eigenvalue problem for the most recent scalar/AO Fock matrices. Add the associated coefficient matrices and orbital energies to the environment.
     * 
     *  @param environment              The environment that acts as a sort of calculation space.
     */
    void execute(Environment& environment) override {

        const auto& F = environment.fock_matrices.back();  // The most recent scalar/AO basis alpha & beta Fock matrix.

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver_alpha {F.alpha().parameters(), environment.S.alpha().parameters()};
        const UTransformationComponent<Scalar>& C_alpha {generalized_eigensolver_alpha.eigenvectors()};

        const auto& orbital_energies_alpha = generalized_eigensolver_alpha.eigenvalues();

        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver_beta {F.beta().parameters(), environment.S.beta().parameters()};
        const UTransformationComponent<Scalar>& C_beta {generalized_eigensolver_beta.eigenvectors()};

        const auto& orbital_energies_beta = generalized_eigensolver_beta.eigenvalues();

        const UTransformation<Scalar>& C {C_alpha, C_beta};
        const SpinResolved<VectorX<Scalar>> mo_energies {orbital_energies_alpha, orbital_energies_beta};

        environment.coefficient_matrices.push_back(C);
        environment.orbital_energies.push_back(mo_energies);
    }
};


}  // namespace GQCP
