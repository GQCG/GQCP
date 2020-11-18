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


#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  The generalized Hartree-Fock stability matrices.
 * 
 *  @tparam _Scalar             The type of scalar that is used for the elements of the stability matrices: real or complex.
 */
template <typename _Scalar>
class GHFStabilityMatrices {
public:
    // The type of scalar that is used for the elements of the stability matrices: real or complex.
    using Scalar = _Scalar;

    // The type of one of the two components
    using Matrix = MatrixX<Scalar>;

private:
    // The A submatrix.
    Matrix A;

    // The B submatrix.
    Matrix B;

public:
    /*
     *  MARK: Constructors
     */

    /*
     *  Construct the object containing all building blocks for the GHF stability matrices.
     * 
     *  @param A        The submatrix A.
     *  @param B        The submatrix B.
     */
    GHFStabilityMatrices(const Matrix& A, const Matrix& B) :
        A {A},
        B {B} {}

public:
    /*
     *  MARK: Accessing the submatrices
     */

    /**
     *  @return A read-only reference to the A submatrix.
     */
    const Matrix& subMatrixA() const { return this->A; }

    /**
     *  @return A writable reference to the A submatrix.
     */
    Matrix& subMatrixA() { return this->A; }

    /**
     *  @return A read-only reference to the B submatrix.
     */
    const Matrix& subMatrixB() const { return this->B; }

    /**
     *  @return A writable reference to the B submatrix.
     */
    Matrix& subMatrixB() { return this->B; }


    /*
     *  MARK: Constructing the stability matrices
     */

    /**
     *  Construct the internal stability matrix of the real GHF method.
     *
     *  @note The internal stability condition of the real GHF method is checked using A+B.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> internal() const { return this->subMatrixA() + this->subMatrixB(); }


    /**
     *  Construct the external stability matrix of the real GHF method.
     *
     *  @note The external stability condition of the real GHF method is checked using A-B.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> external() const { return this->subMatrixA() - this->subMatrixB(); }


    /**
     *  @return the internal stability matrix of the complex GHF method.
     *
     *  @note The internal stability condition of the real GHF method is checked using (A,   B  )
     *                                                                                 (B^*, A^*)
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, MatrixX<complex>> internal() const {

        // Calculate the necessary partial stability matrices.
        const auto A = this->subMatrixA();
        const auto B = this->subMatrixB();

        // Determine the dimensions of the total stability matrix.
        const auto K = A.dimension(0);
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        Matrix H {dim, dim};

        H.topLeftCorner(K, K) = A;
        H.topRightCorner(K, K) = B;
        H.bottomLeftCorner(K, K) = B.conjugate();
        H.bottomRightCorner(K, K) = A.conjugate();

        return H;
    }


    /*
     *  MARK: Checking the stability
     */

    /**
     *  @return a boolean, telling us whether or not the real valued internal stability matrix belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, bool> isInternallyStable() const {

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real valued wavefunction.
        const auto stability_matrix = this->internal();

        // Create an eigensolver to diagonalize the stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {stability_matrix};

        // Calculate the eigenvalues and check whether they are strictly positive or not, we update the stability properties accordingly.
        const auto& eigenvalues = eigensolver.eigenvalues().transpose();

        if (eigenvalues[0] < -1.0e-5) {
            return false;
        } else {
            return true;
        }
    }


    /**
     *  @return a boolean, telling us whether or not the real valued external stability matrix belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, bool> isExternallyStable() const {

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real valued wavefunction.
        const auto stability_matrix = this->external();

        // Create an eigensolver to diagonalize the stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {stability_matrix};

        // Calculate the eigenvalues and check whether they are strictly positive or not, we update the stability properties accordingly.
        const auto& eigenvalues = eigensolver.eigenvalues().transpose();

        if (eigenvalues[0] < -1.0e-5) {
            return false;
        } else {
            return true;
        }
    }


    /**
     *  @return a boolean, telling us whether or not the complex valued internal stability matrix belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, bool> isInternallyStable() const {

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real valued wavefunction.
        const auto stability_matrix = this->internal();

        // Create an eigensolver to diagonalize the stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {stability_matrix};

        // Calculate the eigenvalues and check whether they are strictly positive or not, we update the stability properties accordingly.
        const auto& eigenvalues = eigensolver.eigenvalues().transpose();

        if (eigenvalues[0] < -1.0e-5) {
            return false;
        } else {
            return true;
        }
    }


    /*
     *  MARK: printing the stability properties of these stability matrices 
     */

    /*
     *  Print the description of the stability properties of a real valued GHF wavefunction.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The real valued GHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The real valued GHF wavefunction contains an internal instability." << std::endl;
        }

        // Next check the external stability
        if (this->isExternallyStable() == true) {
            std::cout << "The real valued GHF wavefunction is externally stable." << std::endl;
        } else {
            std::cout << "The real valued GHF wavefunction contains a real->complex external instability." << std::endl;
        }
    }


    /*
     *  Print the description of the stability properties of a complex valued GHF wavefunction.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The complex valued GHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The real complex GHF wavefunction contains an internal instability." << std::endl;
        }
    }
};

}  // namespace GQCP
