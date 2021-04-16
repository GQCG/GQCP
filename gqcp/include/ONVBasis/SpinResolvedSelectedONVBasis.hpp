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


#include "Mathematical/Representation/MatrixRepresentationEvaluationContainer.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"


namespace GQCP {


/**
 *  A spin-resolved basis that is flexible in its compromising (spin-resolved) ONVs.
 */
class SpinResolvedSelectedONVBasis {
private:
    // The number of spin-orbitals (equal for alpha and beta).
    size_t K;

    // The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
    size_t N_alpha;

    // The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
    size_t N_beta;

    // A collection of ONVs that span a 'selected' part of a Fock space.
    std::vector<SpinResolvedONV> onvs;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct an empty spin-resolved selected ONV basis.
     *
     *  @param K            The number of spin-orbitals (equal for alpha and beta).
     *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     */
    SpinResolvedSelectedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta);

    /**
     *  The default constructor.
     */
    SpinResolvedSelectedONVBasis() = default;

    /**
     *  Generate a `SpinResolvedSelectedONVBasis` from a seniority-zero ONV basis.
     *
     *  @param onv_basis        The seniority-zero ONV basis.
     */
    SpinResolvedSelectedONVBasis(const SeniorityZeroONVBasis& onv_basis);

    /**
     *  Generate a `SpinResolvedSelectedONVBasis` from a full spin-resolved ONV basis.
     *
     *  @param onv_basis        The full spin-resolved ONV basis.
     */
    SpinResolvedSelectedONVBasis(const SpinResolvedONVBasis& onv_basis);


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a `SpinResolvedSelectedONVBasis` for a CI singles calculation, using the HF determinant as a reference.
     * 
     *  @param K            The number of spin-orbitals (equal for alpha and beta).
     *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     * 
     *  @returen A CI singles-equivalent `SpinResolvedSelectedONVBasis`.
     */
    static SpinResolvedSelectedONVBasis CIS(const size_t K, const size_t N_alpha, const size_t N_beta);


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of spin-orbitals (equal for alpha and beta).
     */
    size_t numberOfOrbitals() const { return this->K; }

    /**
     *  @return The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
     */
    size_t numberOfAlphaElectrons() const { return this->N_alpha; }

    /**
     *  @return The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
     */
    size_t numberOfBetaElectrons() const { return this->N_beta; }

    /**
     *  @return The dimension of the Fock subspace that is spanned by this selected ONV basis.
     */
    size_t dimension() const { return this->onvs.size(); }


    /*
     *  MARK: Modifying
     */

    /**
     *  Expand this ONV basis with the given spin-resolved ONV.
     * 
     *  @param onv          The ONV that should be included in this ONV basis.
     */
    void expandWith(const SpinResolvedONV& onv);

    /**
     *  Expand this ONV basis with the given spin-resolved ONVs.
     * 
     *  @param onvs         The ONVs that should be included in this ONV basis.
     */
    void expandWith(const std::vector<SpinResolvedONV>& onvs);


    /*
     *  MARK: Accessing
     */

    /**
     *  Access the ONV that corresponds to the given index/address.
     * 
     *  @param index            The address of the ONV.
     * 
     *  @return The ONV that corresponds to the given index/address.
     */
    const SpinResolvedONV& onvWithIndex(const size_t index) const { return this->onvs[index]; }


    /*
     *  MARK: Dense restricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarRSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the dense matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the two-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the dense matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Diagonal restricted operator evaluations
     */

    /**
     *  Calculate the diagonal of the matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f_op             A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the one-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const;

    /**
     *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g_op             A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the two-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g_op) const;

    /**
     *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Sparse restricted operator evaluations
     */

    /**
     *  Calculate the sparse matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the one-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarRSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the sparse matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the two-electron operator.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarRSQTwoElectronOperator<double>& g) const;

    /**
     *  Calculate the sparse matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the Hamiltonian.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const RSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Restricted matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted one-electron operator with the given coefficient vector.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted two-electron operator with the given coefficient vector.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& g, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;


    /*
     *  MARK: Dense unrestricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of an unrestricted one-electron operator in this ONV basis.
     *
     *  @param f                An unrestricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the one-electron operator.
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarUSQOneElectronOperator<double>& f) const;

    /**
     *  Calculate the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Diagonal unrestricted operator evaluations
     */

    /**
     *  Calculate the diagonal of the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Sparse unrestricted operator evaluations
     */

    /**
     *  Calculate the sparse matrix representation of an unrestricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A sparse matrix represention of the Hamiltonian.
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const USQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Unrestricted matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) an unrestricted Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x) const;


    /*
     *  MARK: Operator evaluations - general implementations - containers
     */

    /**
     *  Calculate the matrix representation of an unrestricted one-electron operator in this ONV basis and emplace it in the given container.
     * 
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *
     *  @param f                            An unrestricted one-electron operator expressed in an orthonormal spin-orbital basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix>
    void evaluate(const ScalarUSQOneElectronOperator<double>& f, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        const auto dim = this->dimension();
        const auto& f_a = f.alpha().parameters();
        const auto& f_b = f.beta().parameters();

        for (; !container.isFinished(); container.increment()) {
            SpinResolvedONV onv_I = this->onvWithIndex(container.index);
            SpinUnresolvedONV alpha_I = onv_I.onv(Spin::alpha);
            SpinUnresolvedONV beta_I = onv_I.onv(Spin::beta);

            // Calculate the diagonal elements.
            for (size_t p = 0; p < this->numberOfOrbitals(); p++) {
                if (alpha_I.isOccupied(p)) {
                    container.addRowwise(container.index, f_a(p, p));
                }

                if (beta_I.isOccupied(p)) {
                    container.addRowwise(container.index, f_b(p, p));
                }
            }

            // Calculate the off-diagonal elements, by going over all other ONVs J. (I != J)
            for (size_t J = container.index + 1; J < dim; J++) {

                SpinResolvedONV onv_J = this->onvWithIndex(J);
                SpinUnresolvedONV alpha_J = onv_J.onv(Spin::alpha);
                SpinUnresolvedONV beta_J = onv_J.onv(Spin::beta);

                // 1 excitation in the alpha part, 0 in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other.
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.

                    // Calculate the total sign and emplace the evaluation in the container.
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);
                    const double value = f_a(p, q);

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);
                }

                // 0 excitations in alpha part, 1 in the beta.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other.
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.

                    // Calculate the total sign and emplace the evaluation in the container.
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);
                    const double value = f_b(p, q);

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);
                }
            }  // loop over addresses J > I
        }      // container loop
    }


    /**
     *  Calculate the matrix representation of an unrestricted Hamiltonian in this ONV basis and emplace it in the given container.
     *
     *  @tparam Matrix                      The type of matrix used to store the evaluations.
     *
     *  @param hamiltonian                  An unrestricted Hamiltonian expressed in an orthonormal spin-orbital basis.
     *  @param container                    A specialized container for emplacing evaluations/matrix elements.
     */
    template <typename Matrix>
    void evaluate(const USQHamiltonian<double>& hamiltonian, MatrixRepresentationEvaluationContainer<Matrix>& container) const {

        // Prepare some variables.
        const size_t dim = this->dimension();
        const size_t K = this->numberOfOrbitals();

        const auto& h_a = hamiltonian.core().alpha().parameters();
        const auto& g_aa = hamiltonian.twoElectron().alphaAlpha().parameters();
        const auto& h_b = hamiltonian.core().beta().parameters();
        const auto& g_bb = hamiltonian.twoElectron().betaBeta().parameters();

        // For the mixed two-electron integrals g_ab and g_ba, we can use the following relation: g_ab(pqrs) = g_ba(rspq) and proceed to only work with g_ab.
        const auto& g_ab = hamiltonian.twoElectron().alphaBeta().parameters();

        for (; !container.isFinished(); container.increment()) {  // loop over all addresses (I)
            SpinResolvedONV onv_I = this->onvWithIndex(container.index);
            SpinUnresolvedONV alpha_I = onv_I.onv(Spin::alpha);
            SpinUnresolvedONV beta_I = onv_I.onv(Spin::beta);

            // Calculate the diagonal elements (I=J).
            for (size_t p = 0; p < K; p++) {
                if (alpha_I.isOccupied(p)) {
                    container.addRowwise(container.index, h_a(p, p));
                    for (size_t q = 0; q < K; q++) {

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (alpha_I.isOccupied(q)) {
                                container.addRowwise(container.index, 0.5 * g_aa(p, p, q, q));
                                container.addRowwise(container.index, -0.5 * g_aa(p, q, q, p));
                            }
                        }

                        if (beta_I.isOccupied(q)) {
                            container.addRowwise(container.index, 0.5 * g_ab(p, p, q, q));
                        }
                    }  // loop over q
                }

                if (beta_I.isOccupied(p)) {
                    container.addRowwise(container.index, h_b(p, p));
                    for (size_t q = 0; q < K; q++) {

                        if (p != q) {  // can't create/annihilate the same orbital twice
                            if (beta_I.isOccupied(q)) {
                                container.addRowwise(container.index, 0.5 * g_bb(p, p, q, q));
                                container.addRowwise(container.index, -0.5 * g_bb(p, q, q, p));
                            }
                        }

                        if (alpha_I.isOccupied(q)) {
                            container.addRowwise(container.index, 0.5 * g_ab(q, q, p, p));  // g_ab(pqrs) = g_ba(rspq)
                        }
                    }  // loop over q
                }
            }  // loop over q

            // Calculate the off-diagonal elements, by going over all other ONVs (J>I).
            for (size_t J = container.index + 1; J < dim; J++) {

                SpinResolvedONV onv_J = this->onvWithIndex(J);
                SpinUnresolvedONV alpha_J = onv_J.onv(Spin::alpha);
                SpinUnresolvedONV beta_J = onv_J.onv(Spin::beta);

                // 1 excitation in the alpha part, 0 excitations in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other.
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // We're sure that there is only 1 element in the std::vector<size_t>.

                    // Calculate the total sign and emplace the one-electron contribution in the container.
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);
                    const double value = h_a(p, q);

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);


                    for (size_t r = 0; r < K; r++) {                           // r loops over spatial orbitals
                        if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {                        // can't create or annihilate the same orbital

                                double value = 0.5 * (g_aa(p, q, r, r) - g_aa(r, q, p, r) - g_aa(p, r, r, q) + g_aa(r, r, p, q));

                                container.addColumnwise(J, sign * value);
                                container.addRowwise(J, sign * value);
                            }
                        }

                        if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the upper-level if-branch

                            double value = 0.5 * 2 * g_ab(p, q, r, r);  // g_ab(pqrs) = g_ba(rspq)

                            container.addColumnwise(J, sign * value);
                            container.addRowwise(J, sign * value);
                        }
                    }
                }

                // 0 excitations in the alpha part, 1 excitation in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {


                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);

                    double value = h_b(p, q);

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);

                    for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                        if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {                      // can't create or annihilate the same orbital
                                double value = 0.5 * (g_bb(p, q, r, r) - g_bb(r, q, p, r) - g_bb(p, r, r, q) + g_bb(r, r, p, q));

                                container.addColumnwise(J, sign * value);
                                container.addRowwise(J, sign * value);
                            }
                        }

                        if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                            double value = 0.5 * 2 * g_ab(r, r, p, q);  // g_ab(pqrs) = g_ba(rspq)

                            container.addColumnwise(J, sign * value);
                            container.addRowwise(J, sign * value);
                        }
                    }
                }

                // 1 excititation in the alpha part, 1 excitation in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    size_t r = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t s = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                    double value = 0.5 * 2 * g_ab(p, q, r, s);  // g_ab(pqrs) = g_ba(rspq)

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);
                }

                // 2 excitations in the alpha part, 0 excitations in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 4) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = alpha_I.findDifferentOccupations(alpha_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = alpha_J.findDifferentOccupations(alpha_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];

                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);

                    double value = 0.5 * (g_aa(p, q, r, s) - g_aa(p, s, r, q) - g_aa(r, q, p, s) + g_aa(r, s, p, q));

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);
                }

                // 0 excitations in the alpha part, 2 excitations in the beta part.
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 4)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = beta_I.findDifferentOccupations(beta_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = beta_J.findDifferentOccupations(beta_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];

                    int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);

                    double value = 0.5 * (g_bb(p, q, r, s) - g_bb(p, s, r, q) - g_bb(r, q, p, s) + g_bb(r, s, p, q));

                    container.addColumnwise(J, sign * value);
                    container.addRowwise(J, sign * value);
                }
            }  // loop over addresses J > I
        }      // loop over addresses I
    }
};


}  // namespace GQCP
