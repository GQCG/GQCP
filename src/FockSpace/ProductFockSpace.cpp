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

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


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
{
    this->alpha_couplings = this->fock_space_alpha.calculateOneElectronCouplings();
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
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("ProductFockSpace::calculateDimension(size_t, size_t, size_t): "+ std::string(e.what()));

    }
}



/*
 * PUBLIC METHODS
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


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const OneElectronOperator<double>& one_op,
                                                             bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(one_op, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(one_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (size_t i = 0; i < alpha_evaluation.cols(); i++){
        for (size_t j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const OneElectronOperator<double>& one_op,
                                                                     bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(OneElectronOperator<double>, bool): Not implemented.");
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const TwoElectronOperator<double>& two_op,
                                                             bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(two_op, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(two_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = this->fock_space_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

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
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const TwoElectronOperator<double>& two_op,
                                                                     bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(TwoElectronOperator<double>, bool): Not implemented.");
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param ham_par              Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const HamiltonianParameters<double>& ham_par,
                                                             bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(ham_par, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(ham_par, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, ham_par.get_g());
        const auto& beta_two_electron_intermediate = this->fock_space_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, ham_par.get_g());
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

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
 *  @param ham_par              Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const HamiltonianParameters<double>& ham_par,
                                                                     bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(HamiltonianParameters<double>, bool): Not implemented.");
}


/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const OneElectronOperator<double>& one_op) const {

    auto K = one_op.get_K();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(OneElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);
                diagonal(Ia * dim_beta + Ib) += one_op(p, p);

            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += one_op(p, p);
            }

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const TwoElectronOperator<double>& two_op) const {

    auto K = two_op.get_K();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(TwoElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    OneElectronOperator<double> k = two_op.effectiveOneElectronPartition();

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_alpha.isOccupied(q)) {  // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += two_op(p, p, q, q);
                    }
                }  // q loop
            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op(p, q, q, p);
                    }
                }  // q loop
            }  // e_b loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the Hamiltonian in this Fock space
 *
 *  @param ham_par              Hamiltonian parameters in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const HamiltonianParameters<double>& ham_par) const {
    return this->evaluateOperatorDiagonal(ham_par.get_h()) + this->evaluateOperatorDiagonal(ham_par.get_g());
}



}  // namespace GQCP
