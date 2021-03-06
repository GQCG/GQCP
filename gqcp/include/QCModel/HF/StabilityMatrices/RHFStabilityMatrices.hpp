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
 *  The restricted Hartree-Fock stability matrices.
 * 
 *  @tparam _Scalar             The type of scalar that is used for the elements of the stability matrices: real or complex.
 */
template <typename _Scalar>
class RHFStabilityMatrices {
public:
    // The type of scalar that is used for the elements of the stability matrices: real or complex.
    using Scalar = _Scalar;

    // The type of one of the components
    using Matrix = MatrixX<Scalar>;

private:
    // The singlet A submatrix.
    Matrix singlet_A;

    // The singlet B submatrix.
    Matrix singlet_B;

    // The triplet A submatrix.
    Matrix triplet_A;

    // The triplet B submatrix.
    Matrix triplet_B;

public:
    /*
     *  MARK: Constructors
     */

    /*
     *  Construct the object containing all building blocks for the RHF stability matrices.
     * 
     *  @param singlet_A        The singlet A submatrix.
     *  @param singlet_B        The singlet B submatrix.
     *  @param triplet_A        The triplet A submatrix.
     *  @param triplet_B        The triplet B submatrix.
     */
    RHFStabilityMatrices(const Matrix& singlet_A, const Matrix& singlet_B, const Matrix& triplet_A, const Matrix& triplet_B) :
        singlet_A {singlet_A},
        singlet_B {singlet_B},
        triplet_A {triplet_A},
        triplet_B {triplet_B} {}

public:
    /*
     *  MARK: Accessing the submatrices
     */

    /**
     *  @return A read-only reference to the singlet A submatrix.
     */
    const Matrix& singletA() const { return this->singlet_A; }

    /**
     *  @return A read-only reference to the singlet B submatrix.
     */
    const Matrix& singletB() const { return this->singlet_B; }

    /**
     *  @return A read-only reference to the triplet A submatrix.
     */
    const Matrix& tripletA() const { return this->triplet_A; }

    /**
     *  @return A read-only reference to the triplet B submatrix.
     */
    const Matrix& tripletB() const { return this->triplet_B; }


    /*
     *  MARK: Constructing the stability matrices
     */

    /**
     *  @return The internal stability matrix of the real RHF method.
     *
     *  @note The internal stability condition of the real RHF method is checked using singlet_A + singlet_B.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> internal() const { return this->singletA() + this->singletB(); }


    /**
     *  @return The real->complex external stability matrix of the real RHF method.
     *
     *  @note The real->complex external stability condition of the real RHF method is checked using singlet_A - singlet_B.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> realComplex() const { return this->singletA() - this->singletB(); }


    /**
     *  @return The restricted->unrestricted external stability matrix of the real RHF method.
     *
     *  @note The restricted->unrestricted external stability condition of the real RHF method is checked using triplet_A + triplet_B.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> restrictedUnrestricted() const { return this->tripletA() + this->tripletB(); }


    /**
     *  @return The internal stability matrix of the complex RHF method.
     *
     *  @note The internal stability condition of the complex RHF method is checked using (singlet_A,   singlet_B  )
     *                                                                                    (singlet_B^*, singlet_A^*).
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, MatrixX<complex>> internal() const {

        // Calculate the necessary partial stability matrices.
        const auto& singlet_A = this->singletA();
        const auto& singlet_B = this->singletB();

        // Determine the dimensions of the total stability matrix.
        const auto K = singlet_A.rows();
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        Matrix singlet_H {dim, dim};

        singlet_H.topLeftCorner(K, K) = singlet_A;
        singlet_H.topRightCorner(K, K) = singlet_B;
        singlet_H.bottomLeftCorner(K, K) = singlet_B.conjugate();
        singlet_H.bottomRightCorner(K, K) = singlet_A.conjugate();

        return singlet_H;
    }


    /**
     *  @return The restricted->unrestricted external stability matrix of the complex RHF method.
     *
     *  @note The restricted->unrestricted external stability condition of the complex RHF method is checked using (triplet_A,   triplet_B  )
     *                                                                                                             (triplet_B^*, triplet_A^*).
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, MatrixX<complex>> restrictedUnrestricted() const {

        // Calculate the necessary partial stability matrices.
        const auto& triplet_A = this->tripletA();
        const auto& triplet_B = this->tripletB();

        // Determine the dimensions of the total stability matrix.
        const auto K = triplet_A.rows();
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        Matrix triplet_H {dim, dim};

        triplet_H.topLeftCorner(K, K) = triplet_A;
        triplet_H.topRightCorner(K, K) = triplet_B;
        triplet_H.bottomLeftCorner(K, K) = triplet_B.conjugate();
        triplet_H.bottomRightCorner(K, K) = triplet_A.conjugate();

        return triplet_H;
    }


    /*
     *  MARK: Checking the stability
     */

    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite.
     * 
     *  @return A boolean, telling us if the real or complex valued internal stability matrix belongs to a stable or unstable set of parameters.
     */
    const bool isInternallyStable(const double threshold = -1.0e-5) const {

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real or complex valued wavefunction.
        const auto stability_matrix = this->internal();

        // Check if the stability matrix is positive semi-definite. This indicates stability.
        return stability_matrix.isPositiveSemiDefinite(threshold);
    }


    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite.
     * 
     *  @return A boolean, telling us if the real or complex valued restricted->unrestricted stability matrix belongs to a stable or unstable set of parameters.
     */
    const bool isTripletStable(const double threshold = -1.0e-5) const {

        // The first step is to calculate the correct stability matrix: This method checks the restricted->unrestricted stability of a real or complex valued wavefunction.
        const auto stability_matrix = this->restrictedUnrestricted();

        // Check if the stability matrix is positive semi-definite. This indicates stability.
        return stability_matrix.isPositiveSemiDefinite(threshold);
    }


    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite. 
     *
     *  @return A boolean, telling us if the real->complex stability matrix belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, bool> isComplexConjugateStable(const double threshold = -1.0e-5) const {

        // The first step is to calculate the correct stability matrix: This method checks the real->complex stability of a real valued wavefunction.
        const auto stability_matrix = this->realComplex();

        // Check if the stability matrix is positive semi-definite. This indicates stability.
        return stability_matrix.isPositiveSemiDefinite(threshold);
    }


    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite.
     * 
     *  @return A boolean, telling us whether the parameters are completely externally stable.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, bool> isExternallyStable(const double threshold = -1.0e-5) const {

        // Return whether the parameters are completely externally stable.
        return this->isComplexConjugateStable(threshold) && this->isTripletStable(threshold);
    }


    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite.
     * 
     *  @return A boolean, telling us if the complex valued external stability matrices belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, bool> isExternallyStable(const double threshold = -1.0e-5) const {
        return this->isTripletStable(threshold);
    }


    /*
     *  MARK: printing the stability properties of these stability matrices 
     */

    /*
     *  Print the description of the stability properties of a real valued RHF wavefunction.
     * 
     *  @note   This method runs the stability calculation before printing the results.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The real valued RHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The real valued RHF wavefunction contains an internal instability." << std::endl;
        }

        // The real->complex external stability.
        if (this->isComplexConjugateStable() == true) {
            std::cout << "The real valued RHF wavefunction is stable within the real/complex RHF space." << std::endl;
        } else {
            std::cout << "The real valued RHF wavefunction contains a real->complex instability." << std::endl;
        }

        // The restricted->unrestricted external stability.
        if (this->isTripletStable() == true) {
            std::cout << "The real valued RHF wavefunction is stable within the real RHF/UHF space." << std::endl;
        } else {
            std::cout << "The real valued RHF wavefunction contains a restricted->unrestricted instability." << std::endl;
        }
    }


    /*
     *  Print the description of the stability properties of a complex valued RHF wavefunction.
     * 
     *  @note   This method runs the stability calculation before printing the results.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The complex valued RHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The complex valued RHF wavefunction contains an internal instability." << std::endl;
        }

        // First check the internal stability
        if (this->isExternallyStable() == true) {
            std::cout << "The complex valued RHF wavefunction is stable within the complex RHF/UHF space." << std::endl;
        } else {
            std::cout << "The complex valued RHF wavefunction contains a restricted->unrestricted instability." << std::endl;
        }
    }
};

}  // namespace GQCP
