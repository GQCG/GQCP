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

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Construct an empty spin-resolved selected ONV basis.
 *
 *  @param M            The number of spin-orbitals (equal for alpha and beta).
 *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
 *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta) :
    M {M},
    N_alpha {N_alpha},
    N_beta {N_beta} {}


/**
 *  Generate a `SpinResolvedSelectedONVBasis` from a seniority-zero ONV basis.
 *
 *  @param onv_basis        The seniority-zero ONV basis.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SeniorityZeroONVBasis& onv_basis) :
    SpinResolvedSelectedONVBasis(onv_basis.numberOfSpatialOrbitals(), onv_basis.numberOfElectronPairs(), onv_basis.numberOfElectronPairs()) {

    // Prepare some variables.
    const auto dimension = onv_basis.dimension();
    const auto proxy_onv_basis = onv_basis.proxy();

    // Iterate over the seniority-zero ONV basis and add all ONVs as doubly-occupied ONVs.
    std::vector<SpinResolvedONV> onvs;
    SpinUnresolvedONV onv = proxy_onv_basis.constructONVFromAddress(0);
    for (size_t I = 0; I < dimension; I++) {  // I iterates over all addresses of the doubly-occupied ONVs

        onvs.emplace_back(onv, onv);

        if (I < dimension - 1) {  // prevent the last permutation from occurring
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }
    }

    this->onvs = onvs;
}


/**
 *  Generate a `SpinResolvedSelectedONVBasis` from a full spin-resolved ONV basis.
 *
 *  @param onv_basis        The full spin-resolved ONV basis.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SpinResolvedONVBasis& onv_basis) :
    SpinResolvedSelectedONVBasis(onv_basis.alpha().numberOfOrbitals(), onv_basis.alpha().numberOfElectrons(), onv_basis.beta().numberOfElectrons()) {

    std::vector<SpinResolvedONV> onvs;

    const SpinUnresolvedONVBasis& onv_basis_alpha = onv_basis.alpha();
    const SpinUnresolvedONVBasis& onv_basis_beta = onv_basis.beta();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    SpinUnresolvedONV alpha = onv_basis_alpha.constructONVFromAddress(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        SpinUnresolvedONV beta = onv_basis_beta.constructONVFromAddress(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            onvs.push_back(SpinResolvedONV {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(alpha);
        }
    }
    this->onvs = onvs;
}


/*
 *  MARK: Modifying
 */

/**
 *  Expand this ONV basis with the given spin-resolved ONV.
 * 
 *  @param onv          The ONV that should be included in this ONV basis.
 */
void SpinResolvedSelectedONVBasis::expandWith(const SpinResolvedONV& onv) {

    this->onvs.push_back(onv);
}


/**
 *  Expand this ONV basis with the given spin-resolved ONVs.
 * 
 *  @param onvs         The ONVs that should be included in this ONV basis.
 */
void SpinResolvedSelectedONVBasis::expandWith(const std::vector<SpinResolvedONV>& onvs) {

    for (const auto& onv : onvs) {
        this->expandWith(onv);
    }
}


/*
 *  MARK: Operator evaluations
 */


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
 */
// SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarRSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

//     const auto K = one_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(ScalarRSQOneElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator {this->dimension()};
//     this->evaluate<SquareMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the operator in a dense matrix
//  *
//  *  @param two_op               the two-electron operator to be evaluated in this ONV basis
//  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
//  *
//  *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
//  */
// SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const {

//     const auto K = two_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(ScalarRSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator {this->dimension()};
//     this->evaluate<SquareMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a dense matrix
//  *
//  *  @param sq_hamiltonian               HamiltonianParameters to be evaluated in this ONV basis
//  *  @param diagonal_values              bool to indicate if diagonal values will be calculated
//  *
//  *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
//  */
// SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

//     const auto K = sq_hamiltonian.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(RSQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator {this->dimension()};
//     this->evaluate<SquareMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a dense matrix
//  *
//  *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis
//  *  @param diagonal_values          bool to indicate if diagonal values will be calculated
//  *
//  *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
//  */
// SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, const bool diagonal_values) const {

//     const auto K = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
//     }

//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<SquareMatrix<double>> evaluation_iterator {this->dimension()};
//     this->evaluate<SquareMatrix<double>>(usq_hamiltonian, evaluation_iterator, diagonal_values);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the diagonal of the operator in this ONV basis
//  *
//  *  @param one_op               the one-electron operator to be evaluated in this ONV basis
//  *
//  *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& one_op) const {

//     const auto K = one_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(ScalarRSQTwoElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     const auto& one_op_par = one_op.parameters();

//     // Diagonal contributions
//     VectorX<double> diagonal = VectorX<double>::Zero(dim);

//     for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
//         SpinResolvedONV configuration_I = this->onvWithIndex(I);
//         SpinUnresolvedONV alpha_I = configuration_I.onv(Spin::alpha);
//         SpinUnresolvedONV beta_I = configuration_I.onv(Spin::beta);

//         for (size_t p = 0; p < K; p++) {
//             if (alpha_I.isOccupied(p)) {
//                 diagonal(I) += one_op_par(p, p);
//             }

//             if (beta_I.isOccupied(p)) {
//                 diagonal(I) += one_op_par(p, p);
//             }
//         }  // loop over q

//     }  // alpha address (Ia) loop

//     return diagonal;
// };


// /**
//  *  Evaluate the diagonal of the operator in this ONV basis
//  *
//  *  @param two_op               the two-electron operator to be evaluated in this ONV basis
//  *
//  *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& two_op) const {

//     const auto K = two_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(ScalarRSQTwoElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     const auto& two_op_par = two_op.parameters();

//     // Diagonal contributions
//     VectorX<double> diagonal = VectorX<double>::Zero(dim);

//     for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
//         SpinResolvedONV configuration_I = this->onvWithIndex(I);
//         SpinUnresolvedONV alpha_I = configuration_I.onv(Spin::alpha);
//         SpinUnresolvedONV beta_I = configuration_I.onv(Spin::beta);

//         for (size_t p = 0; p < K; p++) {
//             if (alpha_I.isOccupied(p)) {
//                 for (size_t q = 0; q < K; q++) {

//                     if (p != q) {  // can't create/annihilate the same orbital twice
//                         if (alpha_I.isOccupied(q)) {
//                             diagonal(I) += 0.5 * two_op_par(p, p, q, q);
//                             diagonal(I) -= 0.5 * two_op_par(p, q, q, p);
//                         }
//                     }

//                     if (beta_I.isOccupied(q)) {
//                         diagonal(I) += 0.5 * two_op_par(p, p, q, q);
//                     }
//                 }  // loop over q
//             }

//             if (beta_I.isOccupied(p)) {
//                 for (size_t q = 0; q < K; q++) {

//                     if (p != q) {  // can't create/annihilate the same orbital twice
//                         if (beta_I.isOccupied(q)) {
//                             diagonal(I) += 0.5 * two_op_par(p, p, q, q);
//                             diagonal(I) -= 0.5 * two_op_par(p, q, q, p);
//                         }
//                     }

//                     if (alpha_I.isOccupied(q)) {
//                         diagonal(I) += 0.5 * two_op_par(p, p, q, q);
//                     }
//                 }  // loop over q
//             }
//         }  // loop over q

//     }  // alpha address (Ia) loop

//     return diagonal;
// };


// /**
//  *  Evaluate the diagonal of the Hamiltonian in this ONV basis
//  *
//  *  @param sq_hamiltonian           HamiltonianParameters to be evaluated in this ONV basis
//  *
//  *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const {
//     return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
// };


// /**
//  *  Evaluate the diagonal of the Hamiltonian
//  *
//  *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis
//  *
//  *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {

//     const auto K = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>, bool): Different spinor dimensions of spin components are currently not supported.");
//     }

//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     const auto& h_a = usq_hamiltonian.spinHamiltonian(Spin::alpha).core().parameters();
//     const auto& g_a = usq_hamiltonian.spinHamiltonian(Spin::alpha).twoElectron().parameters();
//     const auto& h_b = usq_hamiltonian.spinHamiltonian(Spin::beta).core().parameters();
//     const auto& g_b = usq_hamiltonian.spinHamiltonian(Spin::beta).twoElectron().parameters();

//     // Only g_ab is stored, for integrals derived from g_ba we reverse the indices as follows : g_ab(pqrs) = g_ba(rspq)
//     const auto& g_ab = usq_hamiltonian.twoElectronMixed().parameters();

//     // Diagonal contributions
//     VectorX<double> diagonal = VectorX<double>::Zero(dim);
//     for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
//         SpinResolvedONV configuration_I = this->onvWithIndex(I);
//         SpinUnresolvedONV alpha_I = configuration_I.onv(Spin::alpha);
//         SpinUnresolvedONV beta_I = configuration_I.onv(Spin::beta);

//         for (size_t p = 0; p < K; p++) {
//             if (alpha_I.isOccupied(p)) {

//                 diagonal(I) += h_a(p, p);

//                 for (size_t q = 0; q < K; q++) {

//                     if (p != q) {  // can't create/annihilate the same orbital twice
//                         if (alpha_I.isOccupied(q)) {
//                             diagonal(I) += 0.5 * g_a(p, p, q, q);
//                             diagonal(I) -= 0.5 * g_a(p, q, q, p);
//                         }
//                     }

//                     if (beta_I.isOccupied(q)) {
//                         diagonal(I) += 0.5 * g_ab(p, p, q, q);
//                     }
//                 }  // loop over q
//             }

//             if (beta_I.isOccupied(p)) {

//                 diagonal(I) += h_b(p, p);

//                 for (size_t q = 0; q < K; q++) {

//                     if (p != q) {  // can't create/annihilate the same orbital twice
//                         if (beta_I.isOccupied(q)) {
//                             diagonal(I) += 0.5 * g_b(p, p, q, q);
//                             diagonal(I) -= 0.5 * g_b(p, q, q, p);
//                         }
//                     }

//                     if (alpha_I.isOccupied(q)) {
//                         diagonal(I) += 0.5 * g_ab(q, q, p, p);
//                     }
//                 }  // loop over q
//             }
//         }  // loop over q

//     }  // alpha address (Ia) loop

//     return diagonal;
// }


// /**
//  *  Evaluate a one electron operator in a matrix vector product
//  *
//  *  @param one_op                       the one electron operator expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in this ONV basis
//  *
//  *  @return the one electron operator's matrix vector product in a vector with the dimensions of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
//     auto K = one_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarRSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator {x, diagonal};
//     this->evaluate<VectorX<double>>(one_op, evaluation_iterator, false);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate a two electron operator in a matrix vector product
//  *
//  *  @param two_op                       the two electron operator expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in this ONV basis
//  *
//  *  @return the two electron operator's matrix vector product in a vector with the dimensions of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
//     auto K = two_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarRSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator {x, diagonal};
//     this->evaluate<VectorX<double>>(two_op, evaluation_iterator, false);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a matrix vector product
//  *
//  *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in this ONV basis
//  *
//  *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {
//     auto K = sq_hamiltonian.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(RSQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator {x, diagonal};
//     this->evaluate<VectorX<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, false);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a matrix vector product
//  *
//  *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in this ONV basis
//  *
//  *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
//  */
// VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     const auto K = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
//     }

//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<VectorX<double>> evaluation_iterator {x, diagonal};
//     this->evaluate<VectorX<double>>(usq_hamiltonian, evaluation_iterator, false);
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the operator in a sparse matrix
//  *
//  *  @param one_op               the one-electron operator to be evaluated in this ONV basis
//  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
//  *
//  *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
//  */
// Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQOneElectronOperator<double>& one_op, const bool diagonal_values) const {

//     const auto K = one_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(ScalarRSQOneElectronOperator<double>, bool): Basis functions of the ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator {this->dimension()};

//     // Estimate the memory that is needed for the evaluation
//     size_t memory = dim * this->M * (this->N_alpha + this->N_beta);
//     if (diagonal_values) {
//         memory += this->dimension();
//     }

//     evaluation_iterator.reserve(memory);
//     this->evaluate<Eigen::SparseMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
//     evaluation_iterator.addToMatrix();
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the operator in a sparse matrix
//  *
//  *  @param two_op               the two-electron operator to be evaluated in this ONV basis
//  *  @param diagonal_values      bool to indicate if diagonal values will be calculated
//  *
//  *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
//  */
// Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const {

//     const auto K = two_op.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(ScalarRSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator {this->dimension()};

//     // Estimate the memory that is needed for the evaluation
//     size_t memory = dim * this->M * this->M * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta);
//     if (diagonal_values) {
//         memory += this->dimension();
//     }

//     evaluation_iterator.reserve(memory);
//     this->evaluate<Eigen::SparseMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
//     evaluation_iterator.addToMatrix();
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a sparse matrix
//  *
//  *  @param sq_hamiltonian               HamiltonianParameters to be evaluated in this ONV basis
//  *  @param diagonal_values              bool to indicate if diagonal values will be calculated
//  *
//  *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
//  */
// Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

//     const auto K = sq_hamiltonian.numberOfOrbitals();
//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(RSQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator {this->dimension()};

//     // Estimate the memory that is needed for the evaluation
//     size_t memory = dim * this->M * this->M * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta) / 4;
//     if (diagonal_values) {
//         memory += this->dimension();
//     }

//     evaluation_iterator.reserve(memory);
//     this->evaluate<Eigen::SparseMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
//     evaluation_iterator.addToMatrix();
//     return evaluation_iterator.evaluation();
// }


// /**
//  *  Evaluate the Hamiltonian in a sparse matrix
//  *
//  *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis
//  *  @param diagonal_values          bool to indicate if diagonal values will be calculated
//  *
//  *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
//  */
// Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const USQHamiltonian<double>& usq_hamiltonian, const bool diagonal_values) const {
//     const auto K = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
//     }

//     if (K != this->M) {
//         throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(USQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
//     }

//     MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> evaluation_iterator {this->dimension()};

//     // Estimate the memory that is needed for the evaluation
//     size_t memory = dim * this->M * this->M * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta);
//     if (diagonal_values) {
//         memory += this->dimension();
//     }

//     evaluation_iterator.reserve(memory);
//     this->evaluate<Eigen::SparseMatrix<double>>(usq_hamiltonian, evaluation_iterator, diagonal_values);
//     evaluation_iterator.addToMatrix();
//     return evaluation_iterator.evaluation();
// }


}  // namespace GQCP
