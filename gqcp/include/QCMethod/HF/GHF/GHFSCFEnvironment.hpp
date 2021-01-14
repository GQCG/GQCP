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


#include "Basis/Transformations/GTransformation.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/aliases.hpp"

#include <Eigen/Dense>

#include <deque>


namespace GQCP {


/**
 *  An algorithm environment that can be used with standard GHF SCF solvers.
 * 
 *  We can basically view it as a compile-time type-safe std::map with all possible information that can be encountered in an GHF SCF algorithm.
 * 
 *  @tparam _Scalar             The scalar type that is used for the coefficient matrix/expansion coefficients: real or complex.
 */
template <typename _Scalar>
class GHFSCFEnvironment {
public:
    using Scalar = _Scalar;


public:
    size_t N;  // The total number of electrons.

    std::deque<Scalar> electronic_energies;

    std::deque<VectorX<Scalar>> orbital_energies;

    ScalarGSQOneElectronOperator<Scalar> S;  // The overlap operator (of both scalar (AO) bases), expressed in spin-blocked notation.

    std::deque<GTransformation<Scalar>> coefficient_matrices;
    std::deque<G1DM<Scalar>> density_matrices;                       // Expressed in the scalar (AO) basis.
    std::deque<ScalarGSQOneElectronOperator<Scalar>> fock_matrices;  // Expressed in the scalar (AO) basis.
    std::deque<VectorX<Scalar>> error_vectors;                       // Expressed in the scalar (AO) basis, used when doing DIIS calculations: the real error matrices should be converted to column-major error vectors for the DIIS algorithm to be used correctly.

    GSQHamiltonian<Scalar> sq_hamiltonian;  // The Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis.


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A constructor that initializes the environment with an initial guess for the coefficient matrix.
     * 
     *  @param N                    The total number of electrons.
     *  @param sq_hamiltonian       The Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis.
     *  @param S                    The overlap operator (of both scalar (AO) bases), expressed in spin-blocked notation.
     *  @param C_initial            The initial coefficient matrix.
     */
    GHFSCFEnvironment(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S, const GTransformation<Scalar>& C_initial) :
        N {N},
        S {S},
        sq_hamiltonian {sq_hamiltonian},
        coefficient_matrices {C_initial} {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Initialize a GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.
     * 
     *  @param N                    The total number of electrons.
     *  @param sq_hamiltonian       The Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis.
     *  @param S                    The overlap operator (of both scalar (AO) bases), expressed in spin-blocked notation.
     * 
     *  @return A GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.
     */
    static GHFSCFEnvironment<Scalar> WithCoreGuess(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S) {

        const auto& H_core = sq_hamiltonian.core().parameters();  // Spin-blocked, in AO basis.

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver {H_core, S.parameters()};
        const GTransformation<Scalar> C_initial {generalized_eigensolver.eigenvectors()};

        return GHFSCFEnvironment<Scalar>(N, sq_hamiltonian, S, C_initial);
    }

    /**
     *  Initialize a GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently adding/subtracting a small complex value from the off-diagonal elements.
     * 
     *  @param N                    The total number of electrons.
     *  @param sq_hamiltonian       The Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis.
     *  @param S                    The overlap operator (of both scalar (AO) bases), expressed in spin-blocked notation.
     * 
     *  @return A GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently adding/subtracting a small complex value from the off-diagonal elements.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, complex>::value, GHFSCFEnvironment<complex>> WithComplexlyTransformedCoreGuess(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S) {

        // Set up the lambda function used to transform the coefficient matrix.
        const auto transformation_function = [](SquareMatrix<complex> C_initial) {
            // Define the complex constant used to transform the initial coefficient matrix.
            const complex x {0, 0.1};

            // Add/subtract the small complex value from the off-diagonal elements of the initial coefficient matrix.
            for (size_t i = 0; i < C_initial.cols(); i++) {
                C_initial(0, i) += x;
            }
            for (size_t j = 0; j < C_initial.rows(); j++) {
                C_initial(j, 0) -= x;
            }

            // Return the updated coefficient matrix.
            return C_initial;
        };

        return GHFSCFEnvironment<complex>::WithTransformedCoreGuess(N, sq_hamiltonian, S, transformation_function);
    }


    /**
     *  Initialize a GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently applying the given unary transformation function.
     * 
     *  @param N                            The total number of electrons.
     *  @param sq_hamiltonian               The Hamiltonian expressed in the scalar (AO) basis, resulting from a quantization using a GSpinorBasis.
     *  @param S                            The overlap operator (of both scalar (AO) bases), expressed in spin-blocked notation.
     *  @param transformation_function      A function that transforms the core guess.
     * 
     *  @return A GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently applying the given unary transformation function.
     */
    static GHFSCFEnvironment<Scalar> WithTransformedCoreGuess(const size_t N, const GSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S, const std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>&)>& transformation_function) {

        const auto& H_core = sq_hamiltonian.core().parameters();  // Spin-blocked, in AO basis.

        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver {H_core, S.parameters()};
        auto C_initial {generalized_eigensolver.eigenvectors()};

        const GTransformation<Scalar> C_initial_complex {transformation_function(C_initial)};

        return GHFSCFEnvironment<Scalar>(N, sq_hamiltonian, S, C_initial_complex);
    }
};


}  // namespace GQCP
