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


#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/complex.hpp"


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


    /*
     *  MARK: Following internal instabilities
     */

    /*
     * Generate the rotation matrix of a real valued UHF wavefunction that will lead to a lower lying minimum.
     *
     * @param occupied_alpha_orbitals       The amount of occupied alpha orbitals in the system.
     * @param occupied_beta_orbitals        The amount of occupied beta orbitals in the system.
     * @param virtual_alpha_orbitals        The amount of virtual alpha orbitals in the system.
     * @param virtual_beta_orbitals         The amount of virtual beta orbitals in the system.
     *
     * @return  The transformation that rotates the solution in the direction of the lowest Hessian eigenvector, towards a global minimum.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, double>::value, UTransformation<double>> instabilityRotationMatrix(const size_t occupied_alpha_orbitals, const size_t occupied_beta_orbitals, const size_t virtual_alpha_orbitals, const size_t virtual_beta_orbitals) const {

        // Calculate the internal stability matrix for the real valued GHF wavefunction.
        const auto H = this->internal();

        // Set up the eigensolver to diagonalize the Hessian/Stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {H};

        // Calculate the lowest eigenvetor of the stability matrix. Split it into its alpha and beta components.
        const auto& eigenvectors = eigensolver.eigenvectors();
        const auto& lowest_eigenvector = eigenvectors.col(0);

        std::vector<double> lowest_eigenvector_alpha;
        std::vector<double> lowest_eigenvector_beta;

        for (size_t i = 0; i < lowest_eigenvector.rows(); i++) {
            if (i < occupied_alpha_orbitals * virtual_alpha_orbitals) {
                lowest_eigenvector_alpha.push_back(lowest_eigenvector[i]);
            } else {
                lowest_eigenvector_beta.push_back(lowest_eigenvector[i]);
            }
        }

        // Create an alpha and a beta kappa matrix.
        Matrix kappa_alpha {occupied_alpha_orbitals, virtual_alpha_orbitals};
        Matrix kappa_beta {occupied_beta_orbitals, virtual_beta_orbitals};

        for (int ra = 0; ra < occupied_alpha_orbitals; ra++) {
            for (int ca = 0; ca < virtual_alpha_orbitals; ca++) {

                // The columns are looped first. Hence, the position within the row is determined by c (+c).
                // The rows are determined by the outer loop. When you start filling a new row, you have to skip all the vector elements used for the previous row (r * cols).
                kappa_alpha(ra, ca) = lowest_eigenvector_alpha[ra * virtual_alpha_orbitals + ca];
            }
        }

        for (int rb = 0; rb < occupied_beta_orbitals; rb++) {
            for (int cb = 0; cb < virtual_beta_orbitals; cb++) {
                kappa_beta(rb, cb) = lowest_eigenvector_beta[rb * virtual_beta_orbitals + cb];
            }
        }

        // Define a rotation matrix kappa (for alpha and beta) of dimension (occupied_sigma+virtual_sigma, occupied_sigma+virtual_sigma) and fill it with the correct elements.
        Matrix Ka {occupied_alpha_orbitals + virtual_alpha_orbitals, occupied_alpha_orbitals + virtual_alpha_orbitals};
        Matrix Kb {occupied_beta_orbitals + virtual_beta_orbitals, occupied_beta_orbitals + virtual_beta_orbitals};

        if (occupied_alpha_orbitals != 0) {
            Ka.topLeftCorner(occupied_alpha_orbitals, occupied_alpha_orbitals) = Matrix::Zero(occupied_alpha_orbitals, occupied_alpha_orbitals);
            Ka.topRightCorner(occupied_alpha_orbitals, virtual_alpha_orbitals) = kappa_alpha;
            Ka.bottomLeftCorner(virtual_alpha_orbitals, occupied_alpha_orbitals) = -1 * (kappa_alpha.transpose().conjugate());
            Ka.bottomRightCorner(virtual_alpha_orbitals, virtual_alpha_orbitals) = Matrix::Zero(virtual_alpha_orbitals, virtual_alpha_orbitals);
        } else {
            Ka = Matrix::Zero(occupied_alpha_orbitals + virtual_alpha_orbitals, occupied_alpha_orbitals + virtual_alpha_orbitals);
            for (size_t i = 0; i < Ka.rows(); i++) {
                Ka(i, i) = 1;
            }
        }
        if (occupied_beta_orbitals != 0) {
            Kb.topLeftCorner(occupied_beta_orbitals, occupied_beta_orbitals) = Matrix::Zero(occupied_beta_orbitals, occupied_beta_orbitals);
            Kb.topRightCorner(occupied_beta_orbitals, virtual_beta_orbitals) = kappa_beta;
            Kb.bottomLeftCorner(virtual_beta_orbitals, occupied_beta_orbitals) = -1 * (kappa_beta.transpose().conjugate());
            Kb.bottomRightCorner(virtual_beta_orbitals, virtual_beta_orbitals) = Matrix::Zero(virtual_beta_orbitals, virtual_beta_orbitals);
        } else {
            Kb = Matrix::Zero(occupied_beta_orbitals + virtual_beta_orbitals, occupied_beta_orbitals + virtual_beta_orbitals);
            for (size_t i = 0; i < Ka.rows(); i++) {
                Kb(i, i) = 1;
            }
        }

        return UTransformation<double> {UTransformationComponent<double> {(-1 * Ka).exp()}, UTransformationComponent<double> {(-1 * Kb).exp()}};
    }


    /*
     * Generate the rotation matrix of a real valued UHF wavefunction that will lead to a lower lying minimum.
     *
     * @param occupied_alpha_orbitals       The amount of occupied alpha orbitals in the system.
     * @param occupied_beta_orbitals        The amount of occupied beta orbitals in the system.
     * @param virtual_alpha_orbitals        The amount of virtual alpha orbitals in the system.
     * @param virtual_beta_orbitals         The amount of virtual beta orbitals in the system.
     *
     * @return  The transformation that rotates the solution in the direction of the lowest Hessian eigenvector, towards a global minimum.
     */
    template <typename S = Scalar>
    enable_if_t<std::is_same<S, complex>::value, UTransformation<complex>> instabilityRotationMatrix(const size_t occupied_alpha_orbitals, const size_t occupied_beta_orbitals, const size_t virtual_alpha_orbitals, const size_t virtual_beta_orbitals) const {

        // Calculate the internal stability matrix for the real valued GHF wavefunction.
        const auto H = this->internal();

        // Set up the eigensolver to diagonalize the Hessian/Stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {H};

        // Calculate the lowest eigenvetor of the stability matrix. Split it into its alpha and beta components.
        const auto& eigenvectors = eigensolver.eigenvectors();
        const auto& lowest_eigenvector = eigenvectors.col(0);

        const auto dim = (occupied_alpha_orbitals * virtual_alpha_orbitals) + (occupied_beta_orbitals * virtual_beta_orbitals);

        std::vector<double> lowest_eigenvector_real;
        std::vector<double> lowest_eigenvector_imaginary;

        for (size_t i = 0; i < dim; i++) {
            lowest_eigenvector_real.push_back(lowest_eigenvector[i].real());
            lowest_eigenvector_imaginary.push_back(lowest_eigenvector[dim + i].real());
        }

        const auto dim_alpha = occupied_alpha_orbitals * virtual_alpha_orbitals;

        std::vector<complex> lowest_eigenvector_alpha;
        std::vector<complex> lowest_eigenvector_beta;

        for (size_t i = 0; i < dim_alpha; i++) {
            const complex x_alpha {lowest_eigenvector_real[i], lowest_eigenvector_imaginary[i]};
            const complex x_beta {lowest_eigenvector_real[dim_alpha + i], lowest_eigenvector_imaginary[dim_alpha + i]};
            lowest_eigenvector_alpha.push_back(x_alpha);
            lowest_eigenvector_beta.push_back(x_beta);
        }

        // Create an alpha and a beta kappa matrix.
        Matrix kappa_alpha {occupied_alpha_orbitals, virtual_alpha_orbitals};
        Matrix kappa_beta {occupied_beta_orbitals, virtual_beta_orbitals};

        for (int ra = 0; ra < occupied_alpha_orbitals; ra++) {
            for (int ca = 0; ca < virtual_alpha_orbitals; ca++) {

                // The columns are looped first. Hence, the position within the row is determined by c (+c).
                // The rows are determined by the outer loop. When you start filling a new row, you have to skip all the vector elements used for the previous row (r * cols).
                kappa_alpha(ra, ca) = lowest_eigenvector_alpha[ra * virtual_alpha_orbitals + ca];
            }
        }

        for (int rb = 0; rb < occupied_beta_orbitals; rb++) {
            for (int cb = 0; cb < virtual_beta_orbitals; cb++) {
                kappa_beta(rb, cb) = lowest_eigenvector_beta[rb * virtual_beta_orbitals + cb];
            }
        }

        // Define a rotation matrix kappa (for alpha and beta) of dimension (occupied_sigma+virtual_sigma, occupied_sigma+virtual_sigma) and fill it with the correct elements.
        Matrix Ka {occupied_alpha_orbitals + virtual_alpha_orbitals, occupied_alpha_orbitals + virtual_alpha_orbitals};
        Matrix Kb {occupied_beta_orbitals + virtual_beta_orbitals, occupied_beta_orbitals + virtual_beta_orbitals};

        if (occupied_alpha_orbitals != 0) {
            Ka.topLeftCorner(occupied_alpha_orbitals, occupied_alpha_orbitals) = Matrix::Zero(occupied_alpha_orbitals, occupied_alpha_orbitals);
            Ka.topRightCorner(occupied_alpha_orbitals, virtual_alpha_orbitals) = kappa_alpha;
            Ka.bottomLeftCorner(virtual_alpha_orbitals, occupied_alpha_orbitals) = -1 * (kappa_alpha.transpose().conjugate());
            Ka.bottomRightCorner(virtual_alpha_orbitals, virtual_alpha_orbitals) = Matrix::Zero(virtual_alpha_orbitals, virtual_alpha_orbitals);
        } else {
            Ka = Matrix::Zero(occupied_alpha_orbitals + virtual_alpha_orbitals, occupied_alpha_orbitals + virtual_alpha_orbitals);
            for (size_t i = 0; i < Ka.rows(); i++) {
                Ka(i, i) = complex(1, 0);
            }
        }
        if (occupied_beta_orbitals != 0) {
            Kb.topLeftCorner(occupied_beta_orbitals, occupied_beta_orbitals) = Matrix::Zero(occupied_beta_orbitals, occupied_beta_orbitals);
            Kb.topRightCorner(occupied_beta_orbitals, virtual_beta_orbitals) = kappa_beta;
            Kb.bottomLeftCorner(virtual_beta_orbitals, occupied_beta_orbitals) = -1 * (kappa_beta.transpose().conjugate());
            Kb.bottomRightCorner(virtual_beta_orbitals, virtual_beta_orbitals) = Matrix::Zero(virtual_beta_orbitals, virtual_beta_orbitals);
        } else {
            Kb = Matrix::Zero(occupied_beta_orbitals + virtual_beta_orbitals, occupied_beta_orbitals + virtual_beta_orbitals);
            for (size_t i = 0; i < Ka.rows(); i++) {
                Kb(i, i) = complex(1, 0);
            }
        }

        return UTransformation<complex> {UTransformationComponent<complex> {(-1 * Ka).exp()}, UTransformationComponent<complex> {(-1 * Kb).exp()}};
    }
};

}  // namespace GQCP
