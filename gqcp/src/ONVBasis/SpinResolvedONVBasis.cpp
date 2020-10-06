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

#include "ONVBasis/SpinResolvedONVBasis.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param M            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
SpinResolvedONVBasis::SpinResolvedONVBasis(const size_t M, const size_t N_alpha, const size_t N_beta) :
    BaseONVBasis(M, SpinResolvedONVBasis::calculateDimension(M, N_alpha, N_beta)),
    onv_basis_alpha {SpinUnresolvedONVBasis(M, N_alpha)},
    onv_basis_beta {SpinUnresolvedONVBasis(M, N_beta)} {

    this->alpha_couplings = this->onv_basis_alpha.calculateOneElectronCouplings();
}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param M            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 *
 *  @return the dimension of the spin-resolved ONV basis
 */
size_t SpinResolvedONVBasis::calculateDimension(const size_t M, const size_t N_alpha, const size_t N_beta) {

    double alpha_dim = SpinUnresolvedONVBasis::calculateDimension(M, N_alpha);
    double beta_dim = SpinUnresolvedONVBasis::calculateDimension(M, N_beta);

    try {
        return boost::numeric::converter<size_t, double>::convert(alpha_dim * beta_dim);
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("SpinResolvedONVBasis::calculateDimension(size_t, size_t, size_t): " + std::string(e.what()));
    }
}


/*
 *  PUBLIC METHODS
 */


/**
 *  Auxiliary method in order to calculate "theta(pq)",
 *  it returns a partition of a two-electron operator as one-electron operator
 *  where A (i,j) = T (p, q, i, j).
 *
 *  @param p            first fixed index of the two-electron operator
 *  @param q            second fixed index of the two-electron operator
 *  @param two_op       the two-electron operator
 *
 *  @return a one-electron operator containing a partition of the two-electron operator
 */
ScalarSQOneElectronOperator<double> SpinResolvedONVBasis::calculateOneElectronPartition(size_t p, size_t q, const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto& two_op_par = two_op.parameters();

    const auto M = two_op.numberOfOrbitals();
    SquareMatrix<double> k_par = SquareMatrix<double>::Zero(M);

    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
            k_par(i, j) += two_op_par(p, q, i, j);
        }
    }

    return ScalarSQOneElectronOperator<double>(k_par);
}


/**
 *  Calculate the compound address of an ONV represented by the two given alpha- and beta-addresses.
 * 
 *  @param I_alpha              the alpha-address
 *  @param I_beta               the beta-address
 * 
 *  @return the compound address of an ONV represented by the two given alpha- and beta-addresses.
 */
size_t SpinResolvedONVBasis::compoundAddress(const size_t I_alpha, const size_t I_beta) const {

    const auto dim_beta = this->onv_basis_beta.dimension();

    return I_alpha * dim_beta + I_beta;
}


/**
 *  @return the dimension of this ONV basis
 */
size_t SpinResolvedONVBasis::dimension() const {

    const auto N_alpha = this->onv_basis_alpha.numberOfElectrons();
    const auto N_beta = this->onv_basis_beta.numberOfElectrons();

    return SpinResolvedONVBasis::calculateDimension(this->M, N_alpha, N_beta);
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->dimension());

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    auto beta_evaluation = onv_basis_beta.evaluateOperatorDense(one_op, diagonal_values);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorDense(one_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (size_t i = 0; i < alpha_evaluation.cols(); i++) {
        for (size_t j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i, j) * ones;
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->dimension());

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    auto beta_evaluation = onv_basis_beta.evaluateOperatorDense(two_op, diagonal_values);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorDense(two_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++) {
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i, j) * ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p < this->M; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2];
        const auto& P = this->calculateOneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = this->onv_basis_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
            }
        }

        for (size_t q = p + 1; q < this->M; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p];
            const auto& P = calculateOneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
                for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->dimension());

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    auto beta_evaluation = onv_basis_beta.evaluateOperatorDense(sq_hamiltonian, diagonal_values);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorDense(sq_hamiltonian, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++) {
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i, j) * ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p < this->M; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2];
        const auto& P = this->calculateOneElectronPartition(p, p, sq_hamiltonian.twoElectron());
        const auto& beta_two_electron_intermediate = this->onv_basis_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
            }
        }

        for (size_t q = p + 1; q < this->M; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p];
            const auto& P = calculateOneElectronPartition(p, q, sq_hamiltonian.twoElectron());
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
                for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param diagonal_values                bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, const bool diagonal_values) const {

    const auto M = usq_hamiltonian.numberOfOrbitals() / 2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->dimension());

    auto const& sq_hamiltonian_alpha = usq_hamiltonian.spinHamiltonian(Spin::alpha);
    auto const& sq_hamiltonian_beta = usq_hamiltonian.spinHamiltonian(Spin::beta);
    auto const& mixed_two_electron_operator = usq_hamiltonian.twoElectronMixed();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    auto beta_evaluation = onv_basis_beta.evaluateOperatorDense(sq_hamiltonian_beta, diagonal_values);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorDense(sq_hamiltonian_alpha, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++) {
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i, j) * ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p < this->M; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2];
        const auto& P = this->calculateOneElectronPartition(p, p, mixed_two_electron_operator);
        const auto& beta_two_electron_intermediate = this->onv_basis_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
            }
        }

        for (size_t q = p + 1; q < this->M; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p];
            const auto& P = calculateOneElectronPartition(p, q, mixed_two_electron_operator);
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
                for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the diagonal of the operator in this spin-resolved ONV basis
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const auto M = one_op.numberOfOrbitals();
    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(ScalarSQOneElectronOperator<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    const auto dim_alpha = onv_basis_alpha.dimension();
    const auto dim_beta = onv_basis_beta.dimension();
    const auto& one_op_par = one_op.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    SpinUnresolvedONV onv_alpha = onv_basis_alpha.constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = onv_basis_beta.constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        onv_basis_beta.transformONVCorrespondingToAddress(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < onv_basis_alpha.numberOfElectrons(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.occupationIndexOf(e_a);
                diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);

            }  // e_a loop

            for (size_t e_b = 0; e_b < onv_basis_beta.numberOfElectrons(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.occupationIndexOf(e_b);
                diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);
            }

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the operator in this spin-resolved ONV basis
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto M = two_op.numberOfOrbitals();
    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(ScalarSQTwoElectronOperator<double>): Basis functions of the SpinUnresolvedONV basis and the operator are incompatible.");
    }

    const auto dim_alpha = onv_basis_alpha.dimension();
    const auto dim_beta = onv_basis_beta.dimension();
    const auto& two_op_par = two_op.parameters();
    const auto k = two_op.effectiveOneElectronPartition().parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    SpinUnresolvedONV onv_alpha = onv_basis_alpha.constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = onv_basis_beta.constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        onv_basis_beta.transformONVCorrespondingToAddress(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < onv_basis_alpha.numberOfElectrons(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.occupationIndexOf(e_a);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < this->M; q++) {  // q loops over SOs
                    if (onv_alpha.isOccupied(q)) {      // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += two_op_par(p, p, q, q);
                    }
                }  // q loop
            }      // e_a loop

            for (size_t e_b = 0; e_b < onv_basis_beta.numberOfElectrons(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.occupationIndexOf(e_b);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < this->M; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {       // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
                    }
                }  // q loop
            }      // e_b loop

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the Hamiltonian in this spin-resolved ONV basis
 *
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
}


/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {

    const auto M = usq_hamiltonian.numberOfOrbitals() / 2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Evaluation environment
    auto const& sq_hamiltonian_alpha = usq_hamiltonian.spinHamiltonian(Spin::alpha);
    auto const& sq_hamiltonian_beta = usq_hamiltonian.spinHamiltonian(Spin::beta);
    auto const& mixed_two_electron_operator = usq_hamiltonian.twoElectronMixed();

    const auto dim_alpha = onv_basis_alpha.dimension();
    const auto dim_beta = onv_basis_beta.dimension();
    auto k_alpha = sq_hamiltonian_alpha.core().parameters();
    auto k_beta = sq_hamiltonian_beta.core().parameters();
    const auto& two_op_par_alpha = sq_hamiltonian_alpha.twoElectron().parameters();
    const auto& two_op_par_beta = sq_hamiltonian_beta.twoElectron().parameters();

    k_alpha = k_alpha + sq_hamiltonian_alpha.twoElectron().effectiveOneElectronPartition().parameters();
    k_beta = k_beta + sq_hamiltonian_beta.twoElectron().effectiveOneElectronPartition().parameters();

    // The two_op_par_mixed variable stored as g_aabb, for integrals derived from g_bbaa we reverse the indices as follows : g_aabb(pqrs) = g_bbaa(rspq)
    const auto& two_op_par_mixed = mixed_two_electron_operator.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    SpinUnresolvedONV onv_alpha = onv_basis_alpha.constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = onv_basis_beta.constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        onv_basis_beta.transformONVCorrespondingToAddress(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < onv_basis_alpha.numberOfElectrons(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.occupationIndexOf(e_a);
                diagonal(Ia * dim_beta + Ib) += k_alpha(p, p);

                for (size_t q = 0; q < this->M; q++) {  // q loops over SOs
                    if (onv_alpha.isOccupied(q)) {      // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += two_op_par_mixed(p, p, q, q);
                    }
                }  // q loop
            }      // e_a loop

            for (size_t e_b = 0; e_b < onv_basis_beta.numberOfElectrons(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.occupationIndexOf(e_b);
                diagonal(Ia * dim_beta + Ib) += k_beta(p, p);

                for (size_t q = 0; q < this->M; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {       // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, q, q, p);
                    }
                }  // q loop
            }      // e_b loop

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto M = one_op.numberOfOrbitals();
    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Environment for evaluations
    SpinUnresolvedONVBasis onv_basis_alpha = this->onvBasisAlpha();
    SpinUnresolvedONVBasis onv_basis_beta = this->onvBasisBeta();

    const auto& alpha_couplings = this->alphaCouplings();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    // Map vector to matrix for vectorized multiplications
    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Spin-resolved evaluation
    auto beta_evaluation = onv_basis_beta.evaluateOperatorSparse(one_op, false);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorSparse(one_op, false);

    // Perform the "matvec"
    matvecmap += xmap * alpha_evaluation + beta_evaluation * xmap;

    return matvec;
}


/**
 *  Evaluate a two electron operator in a matrix vector product
 *
 *  @param two_op                       the two electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
 *
 *  @return the two electron operator's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto M = two_op.numberOfOrbitals();
    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Environment for evaluations
    SpinUnresolvedONVBasis onv_basis_alpha = this->onvBasisAlpha();
    SpinUnresolvedONVBasis onv_basis_beta = this->onvBasisBeta();

    const auto& alpha_couplings = this->alphaCouplings();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Mixed-spin evaluation
    for (size_t p = 0; p < this->M; p++) {

        const auto& P = this->calculateOneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorSparse(P, false);

        // matvec : sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2]);
        for (size_t q = p + 1; q < this->M; q++) {

            const auto& P = this->calculateOneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorSparse(P, true);

            // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p]);
        }
    }

    // Spin-resolved evaluation
    auto beta_evaluation = onv_basis_beta.evaluateOperatorSparse(two_op, false);
    auto alpha_evaluation = onv_basis_alpha.evaluateOperatorSparse(two_op, false);

    matvecmap += beta_evaluation * xmap + xmap * alpha_evaluation;

    return matvec;
}


/**
 *  Evaluate the Hamiltonian in a matrix vector product
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
 *
 *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto M = sq_hamiltonian.numberOfOrbitals();
    if (M != this->M) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Environment for evaluations
    const SpinUnresolvedONVBasis& onv_basis_alpha = this->onvBasisAlpha();
    const SpinUnresolvedONVBasis& onv_basis_beta = this->onvBasisBeta();

    const auto& alpha_couplings = this->alphaCouplings();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Mixed-spin evaluation
    for (size_t p = 0; p < this->M; p++) {

        const auto& P = this->calculateOneElectronPartition(p, p, sq_hamiltonian.twoElectron());
        const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorSparse(P, false);

        // matvec : sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2]);
        for (size_t q = p + 1; q < this->M; q++) {

            const auto& P = this->calculateOneElectronPartition(p, q, sq_hamiltonian.twoElectron());
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorSparse(P, true);

            // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p]);
        }
    }

    // Spin-resolved evaluation
    auto beta_hamiltonian = onv_basis_beta.evaluateOperatorSparse(sq_hamiltonian, false);
    auto alpha_hamiltonian = onv_basis_alpha.evaluateOperatorSparse(sq_hamiltonian, false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}


/**
 *  Evaluate the unrestricted Hamiltonian in a matrix vector product
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param x                              the vector upon which the evaluation acts 
 *  @param diagonal                       the diagonal evaluated in the spin-resolved ONV basis
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    auto M = usq_hamiltonian.numberOfOrbitals() / 2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double> , VectorX<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (M != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and usq_hamiltonian are incompatible.");
    }

    // Environment for evaluations
    const SpinUnresolvedONVBasis& onv_basis_alpha = this->onvBasisAlpha();
    const SpinUnresolvedONVBasis& onv_basis_beta = this->onvBasisBeta();

    const auto& alpha_couplings = this->alphaCouplings();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    for (size_t p = 0; p < this->M; p++) {

        const auto& P = this->calculateOneElectronPartition(p, p, usq_hamiltonian.twoElectronMixed());
        const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorDense(P, false);

        // sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2]);
        for (size_t q = p + 1; q < this->M; q++) {

            const auto& P = this->calculateOneElectronPartition(p, q, usq_hamiltonian.twoElectronMixed());
            const auto& beta_two_electron_intermediate = onv_basis_beta.evaluateOperatorDense(P, true);

            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->M + this->M + 1 - p) / 2 + q - p]);
        }
    }

    auto beta_hamiltonian = onv_basis_beta.evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(Spin::beta), false);
    auto alpha_hamiltonian = onv_basis_alpha.evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(Spin::alpha), false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}


/**
 *  Iterate over all ONVs (implicitly, by resolving in their spin components) in this ONV basis and apply the given callback function.
 * 
 *  @param callback             the function to be applied in every iteration. Its arguments are two pairs of spin-unresolved ONVs and their corresponding addresses, where the first two arguments are related to alpha-spin. The last two arguments are related to beta-spin.
 */
void SpinResolvedONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, const size_t, const SpinUnresolvedONV&, const size_t)>& callback) const {

    const auto dim_alpha = this->onv_basis_alpha.dimension();
    const auto dim_beta = this->onv_basis_beta.dimension();

    SpinUnresolvedONV onv_alpha = this->onv_basis_alpha.constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = this->onv_basis_beta.constructONVFromAddress(0);

    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        onv_basis_beta.transformONVCorrespondingToAddress(onv_beta, 0);  // reset the beta ONV to the one with the first address
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {                       // Ib loops over addresses of beta spin strings

            callback(onv_alpha, Ia, onv_beta, Ib);

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop
}


}  // namespace GQCP
