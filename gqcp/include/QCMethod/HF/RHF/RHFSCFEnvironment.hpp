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


#include "Basis/Transformations/RTransformationMatrix.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <Eigen/Dense>

#include <deque>


namespace GQCP {


/**
 *  An algorithm environment that can be used with standard RHF SCF solvers.
 * 
 *  We can basically view it as a compile-time type-safe std::map with all possible information that can be encountered in an RHF SCF algorithm.
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

    std::deque<VectorX<double>> orbital_energies;

    SquareMatrix<Scalar> S;  // the overlap matrix (of the scalar (AO) basis)

    std::deque<RTransformationMatrix<Scalar>> coefficient_matrices;
    std::deque<Orbital1DM<Scalar>> density_matrices;                 // expressed in the scalar (AO) basis
    std::deque<ScalarRSQOneElectronOperator<Scalar>> fock_matrices;  // expressed in the scalar (AO) basis
    std::deque<VectorX<Scalar>> error_vectors;                       // expressed in the scalar (AO) basis, used when doing DIIS calculations: the real error matrices should be converted to column-major error vectors for the DIIS algorithm to be used correctly

    RSQHamiltonian<Scalar> sq_hamiltonian;  // the Hamiltonian expressed in the scalar (AO) basis


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
    RHFSCFEnvironment(const size_t N, const RSQHamiltonian<Scalar>& sq_hamiltonian, const SquareMatrix<Scalar>& S, const RTransformationMatrix<Scalar>& C_initial) :
        N {N},
        S {S},
        sq_hamiltonian {sq_hamiltonian},
        coefficient_matrices {C_initial} {

        if (this->N % 2 != 0) {  // if the total number of electrons is odd
            throw std::invalid_argument("RHFSCFEnvironment::RHFSCFEnvironment(const size_t, const RSQHamiltonian<Scalar>&, const SquareMatrix<Scalar>&, const RTransformationMatrix<Scalar>&): You have given an odd number of electrons.");
        }
    }


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
    static RHFSCFEnvironment<Scalar> WithCoreGuess(const size_t N, const RSQHamiltonian<Scalar>& sq_hamiltonian, const SquareMatrix<Scalar>& S) {

        const auto& H_core = sq_hamiltonian.core().parameters();  // in AO basis

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver {H_core, S};
        const RTransformationMatrix<Scalar> C_initial = generalized_eigensolver.eigenvectors();

        return RHFSCFEnvironment<Scalar>(N, sq_hamiltonian, S, C_initial);
    }
};


}  // namespace GQCP
