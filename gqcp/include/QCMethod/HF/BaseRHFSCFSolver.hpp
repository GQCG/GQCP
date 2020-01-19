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


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"


namespace GQCP {


/**
 *  A base class for the Plain- and DIIS-RHFSCFSolver that implements common properties and methods.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class BaseRHFSCFSolver {
public:
    using ExpansionScalar = _ExpansionScalar;


protected:
    size_t N;  // the total number of electrons
    SquareMatrix<ExpansionScalar> S;  // the overlap matrix of the spinor basis

    SQHamiltonian<ExpansionScalar> sq_hamiltonian;  // expressed in the scalar orbital basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  A general constructor from its properties
     * 
     *  @param N                                    the total number of electrons
     *  @param S                                    the overlap matrix of the spinor basis
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     */
    BaseRHFSCFSolver(const size_t N, const SquareMatrix<ExpansionScalar>& S, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian) :
        N (N),
        S (S),
        sq_hamiltonian (sq_hamiltonian)
    {}


    /*
     *  STATIC PUBLIC METHODS
     */

    static TransformationMatrix<ExpansionScalar> calculateCoreGuess(const ScalarSQOneElectronOperator<ExpansionScalar>& h_core_op, const ScalarSQOneElectronOperator<ExpansionScalar>& S_op) {

        const auto& H_core = h_core_op.parameters();
        const auto S = S_op.parameters();

        // Return an initial guess for the coefficient matrix by solving the generalized eigenvalue problem for H_core
        using MatrixType = Eigen::Matrix<ExpansionScalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (H_core, S);
        TransformationMatrix<ExpansionScalar> C_initial = generalized_eigensolver.eigenvectors();

        return C_initial;
    }
};


}  // namespace GQCP