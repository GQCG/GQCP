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
 *  The unrestricted Hartree-Fock stability matrices.
 * 
 *  @tparam _Scalar             The type of scalar that is used for the elements of the stability matrices: real or complex.
 */
template <typename _Scalar>
class UHFStabilityMatrices {
public:
    // The type of scalar that is used for the elements of the stability matrices: real or complex.
    using Scalar = _Scalar;

    // The type of one of the components.
    using Matrix = MatrixX<Scalar>;

private:
    // The spin-conserved A' submatrix.
    Matrix spin_conserved_A;

    // The spin-conserved B' submatrix.
    Matrix spin_conserved_B;

    // The spin-unconserved A'' submatrix.
    Matrix spin_unconserved_A;

    // The spin-unconserved B'' submatrix.
    Matrix spin_unconserved_B;

public:
    /*
     *  MARK: Constructors
     */

    /*
     *  Construct the object containing all building blocks for the UHF stability matrices.
     * 
     *  @note The names `Spin-conserved` and `Spin-unconserved`come from the article "Constraints and stability in Hartree-Fock theory" by Seeger, R. and Pople J.A. (https://doi.org/10.1063/1.434318).
     * 
     *  @param spin_conserved_A          The spin-conserved A' submatrix.
     *  @param spin_conserved_B          The spin-conserved B' submatrix.
     *  @param spin_unconserved_A        The spin-unconserved A'' submatrix.
     *  @param spin_unconserved_B        The spin-unconserved B'' submatrix.
     */
    UHFStabilityMatrices(const Matrix& spin_conserved_A, const Matrix& spin_conserved_B, const Matrix& spin_unconserved_A, const Matrix& spin_unconserved_B) :
        spin_conserved_A {spin_conserved_A},
        spin_conserved_B {spin_conserved_B},
        spin_unconserved_A {spin_unconserved_A},
        spin_unconserved_B {spin_unconserved_B} {}

public:
    /*
     *  MARK: Accessing the submatrices
     */

    /**
     *  @return A read-only reference to the spin-conserved A' submatrix.
     */
    const Matrix& spinConservedA() const { return this->spin_conserved_A; }

    /**
     *  @return A read-only reference to the spin-conserved B' submatrix.
     */
    const Matrix& spinConservedB() const { return this->spin_conserved_B; }

    /**
     *  @return A read-only reference to the spin-unconserved A'' submatrix.
     */
    const Matrix& spinUnconservedA() const { return this->spin_unconserved_A; }

    /**
     *  @return A read-only reference to the spin-unconserved B'' submatrix.
     */
    const Matrix& spinUnconservedB() const { return this->spin_unconserved_B; }


    /*
     *  MARK: Constructing the stability matrices
     */

    /**
     *  @return The internal stability matrix of the real UHF method.
     *
     *  @note The internal stability condition of the real UHF method is checked using spin-conserved A' + spin-conserved B'.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> internal() const { return this->spinConservedA() + this->spinConservedB(); }


    /**
     *  @return The real->complex external stability matrix of the real UHF method.
     *
     *  @note The real->complex external stability condition of the real UHF method is checked using spin-conserved A' - spin-conserved B'.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> realComplex() const { return this->spinConservedA() - this->spinConservedB(); }


    /**
     *  @return The unrestricted->generalized external stability matrix of the real UHF method.
     *
     *  @note The unrestricted->generalized external stability condition of the real UHF method is checked using spin-unconserved A'' - spin-unconserved B''.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, MatrixX<double>> unrestrictedGeneralized() const { return this->spinUnconservedA() - this->spinUnconservedB(); }


    /**
     *  @return The internal stability matrix of the complex UHF method.
     *
     *  @note The internal stability condition of the complex UHF method is checked using (spin-conserved A',   spin-conserved B'  )
     *                                                                                    (spin-conserved B'^*, spin-conserved A'^*).
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, MatrixX<complex>> internal() const {

        // Calculate the necessary partial stability matrices.
        const auto& spin_conserved_A = this->spinConservedA();
        const auto& spin_conserved_B = this->spinConservedB();

        // Determine the dimensions of the total stability matrix.
        const auto K = spin_conserved_A.rows();
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        Matrix spin_conserved_H {dim, dim};

        spin_conserved_H.topLeftCorner(K, K) = spin_conserved_A;
        spin_conserved_H.topRightCorner(K, K) = spin_conserved_B;
        spin_conserved_H.bottomLeftCorner(K, K) = spin_conserved_B.conjugate();
        spin_conserved_H.bottomRightCorner(K, K) = spin_conserved_A.conjugate();

        return spin_conserved_H;
    }


    /**
     *  @return The unrestricted->generalized external stability matrix of the complex UHF method.
     *
     *  @note The unrestricted->generalized external stability condition of the complex UHF method is checked using (spin-unconserved A',   spin-unconserved B'  )
     *                                                                                                              (spin-unconserved B'^*, spin-unconserved A'^*).
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, MatrixX<complex>> unrestrictedGeneralized() const {

        // Calculate the necessary partial stability matrices.
        const auto& spin_unconserved_A = this->spinUnconservedA();
        const auto& spin_unconserved_B = this->spinUnconservedB();

        // Determine the dimensions of the total stability matrix.
        const auto K = spin_unconserved_A.rows();
        const auto dim = 2 * K;

        // Create the total stability matrix as specified above in the documentation.
        Matrix spin_unconserved_H {dim, dim};

        spin_unconserved_H.topLeftCorner(K, K) = spin_unconserved_A;
        spin_unconserved_H.topRightCorner(K, K) = spin_unconserved_B;
        spin_unconserved_H.bottomLeftCorner(K, K) = spin_unconserved_B.conjugate();
        spin_unconserved_H.bottomRightCorner(K, K) = spin_unconserved_A.conjugate();

        return spin_unconserved_H;
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
     *  @return A boolean, telling us if the real or complex valued unrestricted->generalized stability matrix belongs to a stable or unstable set of parameters.
     */
    const bool isSpinUnconservedStable(const double threshold = -1.0e-5) const {

        // The first step is to calculate the correct stability matrix: This method checks the unrestricted->generalized stability of a real or complex valued wavefunction.
        const auto stability_matrix = this->unrestrictedGeneralized();

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
     *  @return A boolean, telling us whether the real valued parameters are completely externally stable.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, bool> isExternallyStable(const double threshold = -1.0e-5) const {

        // Return whether the parameters are completely externally stable.
        return this->isComplexConjugateStable(threshold) && this->isSpinUnconservedStable(threshold);
    }


    /**
     *  @param threshold        The threshold used to check if the matrix is positive semi-definite. If the lowest eigenvalue is more negative than the threshold, it is not positive semi-definite.
     * 
     *  @return A boolean, telling us if the complex valued external stability matrices belongs to a stable or unstable set of parameters.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, bool> isExternallyStable(const double threshold = -1.0e-5) const {
        return this->isSpinUnconservedStable(threshold);
    }


    /*
     *  MARK: printing the stability properties of these stability matrices 
     */

    /*
     *  Print the description of the stability properties of a real valued UHF wavefunction.
     * 
     *  @note   This method runs the stability calculation before printing the results.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The real valued UHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The real valued UHF wavefunction contains an internal instability." << std::endl;
        }

        // The real->complex external stability.
        if (this->isComplexConjugateStable() == true) {
            std::cout << "The real valued UHF wavefunction is stable within the real/complex UHF space." << std::endl;
        } else {
            std::cout << "The real valued UHF wavefunction contains a real->complex instability." << std::endl;
        }

        // The restricted->unrestricted external stability.
        if (this->isSpinUnconservedStable() == true) {
            std::cout << "The real valued UHF wavefunction is stable within the real UHF/GHF space." << std::endl;
        } else {
            std::cout << "The real valued UHF wavefunction contains an unrestricted->generalized instability." << std::endl;
        }
    }


    /*
     *  Print the description of the stability properties of a complex valued UHF wavefunction.
     * 
     *  @note   This method runs the stability calculation before printing the results.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, void> printStabilityDescription() const {

        // First check the internal stability
        if (this->isInternallyStable() == true) {
            std::cout << "The complex valued UHF wavefunction is internally stable." << std::endl;
        } else {
            std::cout << "The complex valued UHF wavefunction contains an internal instability." << std::endl;
        }

        // First check the internal stability
        if (this->isExternallyStable() == true) {
            std::cout << "The complex valued UHF wavefunction is stable within the complex UHF/GHF space." << std::endl;
        } else {
            std::cout << "The complex valued UHF wavefunction contains an unrestricted->generalized instability." << std::endl;
        }
    }
};

}  // namespace GQCP
