// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "FockSpace/ProductFockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
ProductFockSpace::ProductFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
        BaseFockSpace(K, ProductFockSpace::calculateDimension(K, N_alpha, N_beta)),
        fock_space_alpha (FockSpace(K, N_alpha)),
        fock_space_beta (FockSpace(K, N_beta))
{}



/*
 *  PRIVATE METHODS
 */

OneElectronOperator<double> ProductFockSpace::oneElectronPartition(size_t p, size_t q, const TwoElectronOperator<double>& two_op) const {
    auto K =  two_op.dimension(0);
    OneElectronOperator<double> k = OneElectronOperator<double>::Zero(K, K);

    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < K; j++) {
            k(i, j) += two_op(p, q, i, j);
        }
    }

    return k;
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 *
 *  @return the dimension of the product Fock space
 */
size_t ProductFockSpace::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    double alpha_dim = FockSpace::calculateDimension(K, N_alpha);
    double beta_dim = FockSpace::calculateDimension(K, N_beta);
    try {
        return boost::numeric::converter<size_t, double>::convert(alpha_dim * beta_dim);
    } catch (boost::numeric::bad_numeric_cast &e) {
        throw std::overflow_error("ProductFockSpace::calculateDimension(size_t, size_t, size_t): "+ std::string(e.what()));

    }
}


/*
 * PUBLIC METHODS
 */

/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::EvaluateOperatorDense(const OneElectronOperator<double>& one_op, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorDense(one_op, false);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorDense(one_op, false);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }
    return total_evaluation;
}

/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::EvaluateOperatorSparse(const OneElectronOperator<double>& one_op, bool diagonal_values) const {

    Eigen::SparseMatrix<double> total_evaluation = Eigen::SparseMatrix<double>(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorSparse(one_op, false);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorSparse(one_op, false);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    Eigen::SparseMatrix<double> ones = Eigen::SparseMatrix<double>(dim_beta, dim_beta);
    ones.setIdentity();
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_evaluation, i); it; ++it) {
            total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*ones;
        }
    }
    return total_evaluation;
}

/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::EvaluateOperatorDense(const TwoElectronOperator<double>& two_op, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorDense(two_op, true);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorDense(two_op, true);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = this->fock_space_beta.EvaluateOperatorDense(P, true);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}

/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::EvaluateOperatorSparse(const TwoElectronOperator<double>& two_op, bool diagonal_values) const {
    Eigen::SparseMatrix<double> total_evaluation = Eigen::SparseMatrix<double>(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorSparse(two_op, false);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorSparse(two_op, false);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    Eigen::SparseMatrix<double> ones = Eigen::SparseMatrix<double>(dim_beta, dim_beta);
    ones.setIdentity();
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_evaluation, i); it; ++it) {
            total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = this->fock_space_beta.EvaluateOperatorSparse(P, true);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}

/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::EvaluateOperatorDense(const HamiltonianParameters<double>& ham_par, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorDense(ham_par, true);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorDense(ham_par, true);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, ham_par.get_g());
        const auto& beta_two_electron_intermediate = this->fock_space_beta.EvaluateOperatorDense(P, true);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, ham_par.get_g());
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}

/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::EvaluateOperatorSparse(const HamiltonianParameters<double>& ham_par, bool diagonal_values) const {

    Eigen::SparseMatrix<double> total_evaluation = Eigen::SparseMatrix<double>(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.EvaluateOperatorSparse(ham_par, false);
    auto alpha_evaluation = fock_space_alpha.EvaluateOperatorSparse(ham_par, false);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    Eigen::SparseMatrix<double> ones = Eigen::SparseMatrix<double>(dim_beta, dim_beta);
    ones.setIdentity();
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_evaluation, i); it; ++it) {
            total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, ham_par.get_g());
        const auto& beta_two_electron_intermediate = this->fock_space_beta.EvaluateOperatorSparse(P, true);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, ham_par.get_g());
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}



}  // namespace GQCP
