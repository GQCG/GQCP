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
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  PRIVATE METHODS
 */

/**
 *  @param onv1     the alpha ONV as a string representation read from right to left
 *  @param onv2     the beta ONV as a string representation read from right to left
 *
 *  @return the configuration that holds both ONVs
 *
 *  IMPORTANT: only works for up to 64 bits!
 */
SpinResolvedONV SpinResolvedSelectedONVBasis::makeConfiguration(const std::string& onv1, const std::string& onv2) const {

    boost::dynamic_bitset<> alpha_transfer (onv1);
    boost::dynamic_bitset<> beta_transfer (onv2);

    if (alpha_transfer.size() != this->K | beta_transfer.size() != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::makeConfiguration(std::string, std::string): Given string representations for ONVs are not compatible with the number of orbitals of this ONV basis");
    }

    if (alpha_transfer.count() != this->N_alpha | beta_transfer.count() != this->N_beta) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::makeConfiguration(std::string, std::string): Given string representations for ONVs are not compatible with the number of orbitals of this ONV basis");
    }

    size_t alpha_s = alpha_transfer.to_ulong();
    size_t beta_s = beta_transfer.to_ulong();

    SpinUnresolvedONV alpha (this->K, this->N_alpha, alpha_s);
    SpinUnresolvedONV beta (this->K, this->N_beta, beta_s);

    return SpinResolvedONV {alpha, beta};
}



/*
 *  CONSTRUCTORS
 */

/**
 *  A constructor with initial ONV basis dimension of 0
 *
 *  @param K            the number of orbitals
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(size_t K, size_t N_alpha, size_t N_beta) :
    BaseONVBasis(K, 0),
    N_alpha (N_alpha),
    N_beta (N_beta)
{}


/**
 *  A constructor that generates the onvs based off the given SpinResolvedONVBasis.
 *
 *  @param fock_space       the SpinResolvedONVBasis from which the onvs should be generated
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SpinResolvedONVBasis& fock_space) :
    SpinResolvedSelectedONVBasis (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{
    std::vector<SpinResolvedONV> onvs;

    const SpinUnresolvedONVBasis& fock_space_alpha = fock_space.get_fock_space_alpha();
    const SpinUnresolvedONVBasis& fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    SpinUnresolvedONV alpha = fock_space_alpha.makeONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        SpinUnresolvedONV beta = fock_space_beta.makeONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            onvs.push_back(SpinResolvedONV {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                fock_space_beta.setNextONV(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
            fock_space_alpha.setNextONV(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->onvs = onvs;
}


/**
 *  A constructor that generates 'selected ONVs' based on the given ONVBasis.
 *
 *  @param onv_basis        the seniority-zero ONV basis from which the 'selected ONVs' should be generated
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SeniorityZeroONVBasis& onv_basis) :
    SpinResolvedSelectedONVBasis (onv_basis.numberOfSpatialOrbitals(), onv_basis.numberOfElectronPairs(), onv_basis.numberOfElectronPairs())
{

    // Prepare some variables.
    const auto dimension = onv_basis.dimension();
    const auto proxy_onv_basis = onv_basis.proxy();

    // Iterate over the seniority-zero ONV basis and add all ONVs as doubly-occupied ONVs.
    std::vector<SpinResolvedONV> onvs;
    SpinUnresolvedONV onv = proxy_onv_basis.makeONV(0);
    for (size_t I = 0; I < dimension; I++) {  // I iterates over all addresses of the doubly-occupied ONVs

        onvs.emplace_back(onv, onv);

        if (I < dimension - 1) {  // prevent the last permutation from occurring
            proxy_onv_basis.setNextONV(onv);
        }
    }

    this->dim = dimension;
    this->onvs = onvs;
}


/**
 *  A constructor that generates the onvs based off the given frozen product spin-resolved ONV basis.
 *
 *  @param fock_space       the SpinUnresolvedONVBasis from which the onvs should be generated
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SpinResolvedFrozenONVBasis& fock_space) :
        SpinResolvedSelectedONVBasis (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{
    std::vector<SpinResolvedONV> onvs;

    const SpinUnresolvedFrozenONVBasis& frozen_fock_space_alpha = fock_space.get_frozen_fock_space_alpha();
    const SpinUnresolvedFrozenONVBasis& frozen_fock_space_beta = fock_space.get_frozen_fock_space_beta();

    auto dim_alpha = frozen_fock_space_alpha.get_dimension();
    auto dim_beta = frozen_fock_space_beta.get_dimension();

    SpinUnresolvedONV alpha = frozen_fock_space_alpha.makeONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        SpinUnresolvedONV beta = frozen_fock_space_beta.makeONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            onvs.push_back(SpinResolvedONV {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                frozen_fock_space_beta.setNextONV(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
            frozen_fock_space_alpha.setNextONV(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->onvs = onvs;
}


SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SpinUnresolvedFrozenONVBasis& fock_space) :
        SpinResolvedSelectedONVBasis (fock_space.get_K(), fock_space.get_N(), fock_space.get_N())
{
    std::vector<SpinResolvedONV> onvs;

    auto dim = fock_space.get_dimension();

    // Iterate over the spin-resolved ONV basis and add all onvs as doubly occupied onvs
    SpinUnresolvedONV onv = fock_space.makeONV(0);
    for (size_t I = 0; I < dim; I++) {

        onvs.push_back(SpinResolvedONV {onv, onv});

        if (I < dim - 1) {  // prevent the last permutation from occurring
            fock_space.setNextONV(onv);
        }
    }

    this->dim = dim;
    this->onvs = onvs;
}

/*
 *  PUBLIC METHODS
 */

/**
 *  Make a configuration (see makeConfiguration()) and add it to this spin-resolved ONV basis
 *
 *  @param onv1     the alpha ONV as a string representation read from right to left
 *  @param onv2     the beta ONV as a string representation read from right to left
 */
void SpinResolvedSelectedONVBasis::addConfiguration(const std::string& onv1, const std::string& onv2) {

    this->dim++;

    SpinResolvedONV configuration = makeConfiguration(onv1, onv2);
    onvs.push_back(configuration);
}


/**
 *  Make onvs (see makeConfiguration()) and add them to this ONV basis
 *
 *  @param onv1s     the alpha ONVs as string representations read from right to left
 *  @param onv2s     the beta ONVs as string representations read from right to left
 */
void SpinResolvedSelectedONVBasis::addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s){

    if (onv1s.size() != onv2s.size()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::addConfiguration(const std::string&, const std::string&): Size of both ONV entry vectors do not match");
    }

    for (size_t i = 0; i < onv1s.size(); i++) {
        this->addConfiguration(onv1s[i], onv2s[i]);
    }
}



/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    const auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(ScalarSQOneElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<SquareMatrix<double>> evaluation_iterator (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    const auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(ScalarSQOneElectronOperator<double>, bool): Basis functions of the ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<Eigen::SparseMatrix<double>> evaluation_iterator (this->dim);

    // Estimate the memory that is needed for the evaluation
    size_t memory = dim * this->K * (this->N_alpha + this->N_beta);
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(one_op, evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    const auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(ScalarSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<SquareMatrix<double>> evaluation_iterator (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in this ONV basis
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    const auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(ScalarSQTwoElectronOperator<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<Eigen::SparseMatrix<double>> evaluation_iterator (this->dim);

    // Estimate the memory that is needed for the evaluation
    size_t memory = dim * this->K * this->K * (this->N_alpha + this->N_beta)*(this->N_alpha + this->N_beta);
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(two_op, evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               HamiltonianParameters to be evaluated in this ONV basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {

    const auto K = sq_hamiltonian.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(SQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<SquareMatrix<double>> evaluation_iterator (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param sq_hamiltonian               HamiltonianParameters to be evaluated in this ONV basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {

    const auto K = sq_hamiltonian.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(SQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<Eigen::SparseMatrix<double>> evaluation_iterator (this->dim);

    // Estimate the memory that is needed for the evaluation
    size_t memory = dim * this->K * this->K * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta)/4;
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the diagonal of the operator in this ONV basis
 *
 *  @param one_op               the one-electron operator to be evaluated in this ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(ScalarSQTwoElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    const auto& one_op_par = one_op.parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        SpinResolvedONV configuration_I = this->get_configuration(I);
        SpinUnresolvedONV alpha_I = configuration_I.alphaONV();
        SpinUnresolvedONV beta_I = configuration_I.betaONV();

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                diagonal(I) += one_op_par(p,p);
            }

            if (beta_I.isOccupied(p)) {
                diagonal(I) += one_op_par(p,p);
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;
};

/**
 *  Evaluate the diagonal of the operator in this ONV basis
 *
 *  @param two_op               the two-electron operator to be evaluated in this ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(ScalarSQTwoElectronOperator<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    const auto& two_op_par = two_op.parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        SpinResolvedONV configuration_I = this->get_configuration(I);
        SpinUnresolvedONV alpha_I = configuration_I.alphaONV();
        SpinUnresolvedONV beta_I = configuration_I.betaONV();

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * two_op_par(p,p,q,q);
                            diagonal(I) -= 0.5 * two_op_par(p,q,q,p);
                        }
                    }

                    if (beta_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * two_op_par(p,p,q,q);
                    }
                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * two_op_par(p,p,q,q);
                            diagonal(I) -= 0.5 * two_op_par(p,q,q,p);
                        }
                    }

                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * two_op_par(p,p,q,q);
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;
};


/**
 *  Evaluate the diagonal of the Hamiltonian in this ONV basis
 *
 *  @param sq_hamiltonian           HamiltonianParameters to be evaluated in this ONV basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
};


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<VectorX<double>> evaluation_iterator (x, diagonal);
    this->EvaluateOperator<VectorX<double>>(one_op, evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate a two electron operator in a matrix vector product
 *
 *  @param two_op                       the two electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the two electron operator's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<VectorX<double>> evaluation_iterator (x, diagonal);
    this->EvaluateOperator<VectorX<double>>(two_op, evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a matrix vector product
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = sq_hamiltonian.dimension();
    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<VectorX<double>> evaluation_iterator (x, diagonal);
    this->EvaluateOperator<VectorX<double>>(sq_hamiltonian.core(), sq_hamiltonian.twoElectron(), evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}



/*
 *  UNRESTRICTED
 */ 

/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param diagonal_values          bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
 */
SquareMatrix<double>  SpinResolvedSelectedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const {
   
    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(USQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<SquareMatrix<double>> evaluation_iterator (this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(usq_hamiltonian, evaluation_iterator, diagonal_values);
    return evaluation_iterator.evaluation();

}


/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {

    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>, bool): Different spinor dimensions of spin components are currently not supported.");
    }

    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    const auto& h_a = usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA).core().parameters();
    const auto& g_a = usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA).twoElectron().parameters();
    const auto& h_b = usq_hamiltonian.spinHamiltonian(SpinComponent::BETA).core().parameters();
    const auto& g_b = usq_hamiltonian.spinHamiltonian(SpinComponent::BETA).twoElectron().parameters();

    // Only g_ab is stored, for integrals derived from g_ba we reverse the indices as follows : g_ab(pqrs) = g_ba(rspq)
    const auto& g_ab = usq_hamiltonian.twoElectronMixed().parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);
    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        SpinResolvedONV configuration_I = this->get_configuration(I);
        SpinUnresolvedONV alpha_I = configuration_I.alphaONV();
        SpinUnresolvedONV beta_I = configuration_I.betaONV();

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {

                diagonal(I) += h_a(p,p);

                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g_a(p,p,q,q);
                            diagonal(I) -= 0.5 * g_a(p,q,q,p);
                        }
                    }

                    if (beta_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g_ab(p,p,q,q);
                    }
                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {
    
                diagonal(I) += h_b(p,p);

                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g_b(p,p,q,q);
                            diagonal(I) -= 0.5 * g_b(p,q,q,p);
                        }
                    }

                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g_ab(q,q,p,p);
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the Hamiltonian in a matrix vector product
 *
 *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in this ONV basis
 *
 *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<VectorX<double>> evaluation_iterator (x, diagonal);
    this->EvaluateOperator<VectorX<double>>(usq_hamiltonian, evaluation_iterator, false);
    return evaluation_iterator.evaluation();
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param diagonal_values          bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const {
    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->K) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(USQHamiltonian<double>, bool): Basis functions of this ONV basis and the operator are incompatible.");
    }

    EvaluationIterator<Eigen::SparseMatrix<double>> evaluation_iterator (this->dim);

    // Estimate the memory that is needed for the evaluation
    size_t memory = dim * this->K * this->K * (this->N_alpha + this->N_beta)*(this->N_alpha + this->N_beta);
    if (diagonal_values) {
        memory += this->dim;
    }

    evaluation_iterator.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(usq_hamiltonian, evaluation_iterator, diagonal_values);
    evaluation_iterator.addToMatrix();
    return evaluation_iterator.evaluation();
}


}  // namespace GQCP
