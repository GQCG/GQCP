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
#include "FockSpace/SelectedFockSpace.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


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
Configuration SelectedFockSpace::makeConfiguration(const std::string& onv1, const std::string& onv2) const {

    boost::dynamic_bitset<> alpha_transfer (onv1);
    boost::dynamic_bitset<> beta_transfer (onv2);

    if (alpha_transfer.size() != this->K | beta_transfer.size() != this->K) {
        throw std::invalid_argument("SelectedFockSpace::makeConfiguration(std::string, std::string): Given string representations for ONVs are not compatible with the number of orbitals of the Fock space");
    }

    if (alpha_transfer.count() != this->N_alpha | beta_transfer.count() != this->N_beta) {
        throw std::invalid_argument("SelectedFockSpace::makeConfiguration(std::string, std::string): Given string representations for ONVs are not compatible with the number of orbitals of the Fock space");
    }

    size_t alpha_s = alpha_transfer.to_ulong();
    size_t beta_s = beta_transfer.to_ulong();

    ONV alpha (this->K, this->N_alpha, alpha_s);
    ONV beta (this->K, this->N_beta, beta_s);

    return Configuration {alpha, beta};
}



/*
 *  CONSTRUCTORS
 */

/**
 *  A constructor with initial Fock space dimension of 0
 *
 *  @param K            the number of orbitals
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
SelectedFockSpace::SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
    BaseFockSpace(K, 0),
    N_alpha (N_alpha),
    N_beta (N_beta)
{}


/**
 *  A constructor that generates the configurations based off the given ProductFockSpace.
 *
 *  @param fock_space       the ProductFockSpace from which the configurations should be generated
 */
SelectedFockSpace::SelectedFockSpace(const ProductFockSpace& fock_space) :
    SelectedFockSpace (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{

    std::vector<Configuration> configurations;

    const FockSpace& fock_space_alpha = fock_space.get_fock_space_alpha();
    const FockSpace& fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    ONV alpha = fock_space_alpha.makeONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        ONV beta = fock_space_beta.makeONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            configurations.push_back(Configuration {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation to occur
                fock_space_beta.setNextONV(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            fock_space_alpha.setNextONV(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->configurations = configurations;

}


/**
 *  A constructor that generates the configurations based off the given FockSpace.
 *
 *  @param fock_space       the FockSpace from which the configurations should be generated
 */
SelectedFockSpace::SelectedFockSpace(const FockSpace& fock_space) :
        SelectedFockSpace (fock_space.get_K(), fock_space.get_N(), fock_space.get_N())
{
    std::vector<Configuration> configurations;

    auto dim = fock_space.get_dimension();

    // Iterate over the Fock space and add all onvs as doubly occupied configurations
    ONV onv = fock_space.makeONV(0);
    for (size_t I = 0; I < dim; I++) {

        configurations.push_back(Configuration {onv, onv});

        if (I < dim - 1) {  // prevent the last permutation to occur
            fock_space.setNextONV(onv);
        }
    }

    this->dim = dim;
    this->configurations = configurations;
}


/**
 *  A constructor that generates the configurations based off the given frozen product Fock space.
 *
 *  @param fock_space       the FockSpace from which the configurations should be generated
 */
SelectedFockSpace::SelectedFockSpace(const FrozenProductFockSpace& fock_space) :
        SelectedFockSpace (fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta())
{
    std::vector<Configuration> configurations;

    const FrozenFockSpace& frozen_fock_space_alpha = fock_space.get_frozen_fock_space_alpha();
    const FrozenFockSpace& frozen_fock_space_beta = fock_space.get_frozen_fock_space_beta();

    auto dim_alpha = frozen_fock_space_alpha.get_dimension();
    auto dim_beta = frozen_fock_space_beta.get_dimension();

    ONV alpha = frozen_fock_space_alpha.makeONV(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        ONV beta = frozen_fock_space_beta.makeONV(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            configurations.push_back(Configuration {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation to occur
                frozen_fock_space_beta.setNextONV(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            frozen_fock_space_alpha.setNextONV(alpha);
        }
    }
    this->dim = fock_space.get_dimension();
    this->configurations = configurations;
}


SelectedFockSpace::SelectedFockSpace(const FrozenFockSpace& fock_space) :
        SelectedFockSpace (fock_space.get_K(), fock_space.get_N(), fock_space.get_N())
{
    std::vector<Configuration> configurations;

    auto dim = fock_space.get_dimension();

    // Iterate over the Fock space and add all onvs as doubly occupied configurations
    ONV onv = fock_space.makeONV(0);
    for (size_t I = 0; I < dim; I++) {

        configurations.push_back(Configuration {onv, onv});

        if (I < dim - 1) {  // prevent the last permutation to occur
            fock_space.setNextONV(onv);
        }
    }

    this->dim = dim;
    this->configurations = configurations;
}

/*
 *  PUBLIC METHODS
 */

/**
 *  Make a configuration (see makeConfiguration()) and add it to this Fock space
 *
 *  @param onv1     the alpha ONV as a string representation read from right to left
 *  @param onv2     the beta ONV as a string representation read from right to left
 */
void SelectedFockSpace::addConfiguration(const std::string& onv1, const std::string& onv2) {

    this->dim++;

    Configuration configuration = makeConfiguration(onv1, onv2);
    configurations.push_back(configuration);
}


/**
 *  Make configurations (see makeConfiguration()) and add them to the Fock space
 *
 *  @param onv1s     the alpha ONVs as string representations read from right to left
 *  @param onv2s     the beta ONVs as string representations read from right to left
 */
void SelectedFockSpace::addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s){

    if (onv1s.size() != onv2s.size()) {
        throw std::invalid_argument("SelectedFockSpace::addConfiguration(const std::string&, const std::string&): Size of both ONV entry vectors do not match");
    }

    for (size_t i = 0; i < onv1s.size(); i++) {
        this->addConfiguration(onv1s[i], onv2s[i]);
    }
}



/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> SelectedFockSpace::evaluateOperatorDense(const OneElectronOperator<double>& one_op,
                                                              bool diagonal_values) const {
    EvaluationContainer<SquareMatrix<double>> container(this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(one_op, container, diagonal_values);
    return container.get_container();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> SelectedFockSpace::evaluateOperatorSparse(const OneElectronOperator<double>& one_op,
                                                                      bool diagonal_values) const {
    EvaluationContainer<Eigen::SparseMatrix<double>> container(this->dim);

    // Estimate
    size_t memory = dim * this->K * (this->N_alpha + this->N_beta);
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(one_op, container, diagonal_values);
    container.addToMatrix();
    return container.get_container();
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> SelectedFockSpace::evaluateOperatorDense(const TwoElectronOperator<double>& two_op,
                                                              bool diagonal_values) const {
    EvaluationContainer<SquareMatrix<double>> container(this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(two_op, container, diagonal_values);
    return container.get_container();
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> SelectedFockSpace::evaluateOperatorSparse(const TwoElectronOperator<double>& two_op,
                                                                      bool diagonal_values) const {
    EvaluationContainer<Eigen::SparseMatrix<double>> container(this->dim);

    // Estimate
    size_t memory = dim * this->K * this->K * (this->N_alpha + this->N_beta)*(this->N_alpha + this->N_beta);
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(two_op, container, diagonal_values);
    container.addToMatrix();
    return container.get_container();
}

/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> SelectedFockSpace::evaluateOperatorDense(const HamiltonianParameters<double>& ham_par,
                                                              bool diagonal_values) const {
    EvaluationContainer<SquareMatrix<double>> container(this->dim);
    this->EvaluateOperator<SquareMatrix<double>>(ham_par.get_h(), ham_par.get_g(), container, diagonal_values);
    return container.get_container();
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> SelectedFockSpace::evaluateOperatorSparse(const HamiltonianParameters<double>& ham_par,
                                                                      bool diagonal_values) const {
    EvaluationContainer<Eigen::SparseMatrix<double>> container(this->dim);

    // Estimate
    size_t memory = dim * this->K * this->K * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta)/4;
    if (diagonal_values) {
        memory += this->dim;
    }

    container.reserve(memory);
    this->EvaluateOperator<Eigen::SparseMatrix<double>>(ham_par.get_h(), ham_par.get_g(), container, diagonal_values);

    container.addToMatrix();

    return container.get_container();
}



/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param one_op               the one-electron operator to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> SelectedFockSpace::evaluateOperatorDiagonal(const OneElectronOperator<double>& one_op) const {

    auto K = one_op.cols();
    if (K != this->K) {
        throw std::invalid_argument("SelectedFockSpace::evaluateOperatorDiagonal(TwoElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        Configuration configuration_I = this->get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                diagonal(I) += one_op(p,p);

            }

            if (beta_I.isOccupied(p)) {
                diagonal(I) += one_op(p,p);
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;
};

/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param two_op               the two-electron operator to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> SelectedFockSpace::evaluateOperatorDiagonal(const TwoElectronOperator<double>& two_op) const {

    auto K = two_op.dimension(0);
    if (K != this->K) {
        throw std::invalid_argument("SelectedFockSpace::evaluateOperatorDiagonal(TwoElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        Configuration configuration_I = this->get_configuration(I);
        ONV alpha_I = configuration_I.onv_alpha;
        ONV beta_I = configuration_I.onv_beta;

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * two_op(p,p,q,q);
                            diagonal(I) -= 0.5 * two_op(p,q,q,p);
                        }
                    }

                    if (beta_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * two_op(p,p,q,q);
                    }
                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * two_op(p,p,q,q);
                            diagonal(I) -= 0.5 * two_op(p,q,q,p);
                        }
                    }

                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * two_op(p,p,q,q);
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;

};

/**
 *  Evaluate the diagonal of the Hamiltonian in this Fock space
 *
 *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> SelectedFockSpace::evaluateOperatorDiagonal(const HamiltonianParameters<double>& ham_par) const {
    return this->evaluateOperatorDiagonal(ham_par.get_h()) + this->evaluateOperatorDiagonal(ham_par.get_g());
};



}  // namespace GQCP
