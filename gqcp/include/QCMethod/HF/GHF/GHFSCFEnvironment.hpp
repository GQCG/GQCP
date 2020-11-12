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


#include "Basis/Transformations/GTransformationMatrix.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <Eigen/Dense>

#include <deque>


namespace GQCP {


/**
 *  An algorithm environment that can be used with standard GHF SCF solvers.
 * 
 *  We can basically view it as a compile-time type-safe std::map with all possible information that can be encountered in an GHF SCF algorithm.
 * 
 *  @tparam _Scalar             the scalar type that is used for the coefficient matrix/expansion coefficients
 */
template <typename _Scalar>
class GHFSCFEnvironment {
public:
    using Scalar = _Scalar;


public:
    size_t N;  // the total number of electrons

    std::deque<double> electronic_energies;

    std::deque<VectorX<double>> orbital_energies;

    SquareMatrix<Scalar> S;  // the overlap matrix (of both scalar (AO) bases), expressed in spin-blocked notation

    std::deque<GTransformationMatrix<Scalar>> coefficient_matrices;
    std::deque<G1DM<Scalar>> density_matrices;               // expressed in the scalar (AO) basis
    std::deque<ScalarGSQOneElectronOperator<Scalar>> fock_matrices;  // expressed in the scalar (AO) basis
    std::deque<VectorX<Scalar>> error_vectors;               // expressed in the scalar (AO) basis, used when doing DIIS calculations: the real error matrices should be converted to column-major error vectors for the DIIS algorithm to be used correctly

    GSQHamiltonian<Scalar> sq_hamiltonian;  // the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A constructor that initializes the environment with an initial guess for the coefficient matrix
     * 
     *  @param N                    the total number of electrons
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis
     *  @param S                    the overlap matrix (of both scalar (AO) bases), expressed in spin-blocked notation
     *  @param C_initial            the initial coefficient matrix
     */
    GHFSCFEnvironment(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const SquareMatrix<Scalar>& S, const GTransformationMatrix<Scalar>& C_initial) :
        N {N},
        S {S},
        sq_hamiltonian {sq_hamiltonian},
        coefficient_matrices {C_initial} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Initialize an GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.
     * 
     *  @param N                    the total number of electrons
     *  @param sq_hamiltonian       the Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis
     *  @param S                    the overlap matrix (of both scalar (AO) bases), expressed in spin-blocked notation
     */
    static GHFSCFEnvironment<Scalar> WithCoreGuess(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const SquareMatrix<Scalar>& S) {

        const auto& H_core = sq_hamiltonian.core().parameters();  // spin-blocked, in AO basis

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver {H_core, S};
        const GTransformationMatrix<Scalar> C_initial = generalized_eigensolver.eigenvectors();

        return GHFSCFEnvironment<Scalar>(N, sq_hamiltonian, S, C_initial);
    }
};


}  // namespace GQCP
