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
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/OneRDM.hpp"

#include <Eigen/Dense>


namespace GQCP {


/**
 *  An algorithm environment that can be used with standard RHF SCF solvers.
 * 
 *  @tparam _Scalar             the scalar type that is used for the coefficient matrix/expansion coefficients
 */
template <typename _Scalar>
class RHFSCFEnvironment {
public:
    using Scalar = _Scalar;


public:
    size_t N;  // the total number of electrons

    std::deque<double> electronic_energies;

    QCMatrix<Scalar> S;  // the overlap matrix (of the scalar (AO) basis)
    std::deque<TransformationMatrix<Scalar>> coefficient_matrices;
    std::deque<OneRDM<Scalar>> density_matrices;  // expressed in the scalar (AO) basis
    std::deque<QCMatrix<Scalar>> fock_matrices;  // expressed in the scalar (AO) basis

    SQHamiltonian<Scalar> sq_hamiltonian;  // the Hamiltonian expressed in the scalar (AO) basis


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  A constructor that initializes the environment with an initial guess for the coefficient matrix
     * 
     *  @param N                    the total number of electrons
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis
     *  @param S                    the overlap matrix (of the scalar (AO) basis)
     *  @param C_initial            the initial coefficient matrix
     */
    RHFSCFEnvironment(const size_t N, const SQHamiltonian<Scalar>& sq_hamiltonian, const QCMatrix<Scalar>& S, const TransformationMatrix<Scalar>& C_initial) :
        N (N),
        S (S),
        sq_hamiltonian (sq_hamiltonian),
        coefficient_matrices(1, C_initial)
    {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.
     * 
     *  @param N                    the total number of electrons
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis
     *  @param S                    the overlap matrix (of the scalar (AO) basis)
     */
    RHFSCFEnvironment<Scalar> WithCoreGuess(const size_t N, const SQHamiltonian<Scalar>& sq_hamiltonian, const QCMatrix<Scalar>& S) {

        const auto& H_core = sq_hamiltonian.core().parameters();

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (H_core, S);
        const TransformationMatrix<Scalar> C_initial = generalized_eigensolver.eigenvectors();

        return RHFSCFEnvironment<Scalar>(N, sq_hamiltonian, S, C_initial);
    }
};


}  // namespace GQCP
