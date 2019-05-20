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
#include "HamiltonianBuilder/FCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full alpha and beta product Fock space
 */
FCI::FCI(const ProductFockSpace& fock_space) :
        HamiltonianBuilder(),
        fock_space (fock_space)
{
    FockSpace alpha_fock_space = fock_space.get_fock_space_alpha();
    this->alpha_couplings = alpha_fock_space.calculateOneElectronCouplings();
}

OneElectronOperator<double> FCI::oneElectronPartition(size_t p, size_t q, const TwoElectronOperator<double>& two_op) const {
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
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the FCI Hamiltonian matrix
 */
SquareMatrix<double> FCI::constructHamiltonian(const HamiltonianParameters<double>& hamiltonian_parameters) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("FCI::constructHamiltonian(const HamiltonianParameters<double>&): Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    SquareMatrix<double> total_hamiltonian = SquareMatrix<double>::Zero(this->fock_space.get_dimension(), this->fock_space.get_dimension());
    
    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_hamiltonian = fock_space_beta.EvaluateOperatorSparse(hamiltonian_parameters, false);
    auto alpha_hamiltonian = fock_space_alpha.EvaluateOperatorSparse(hamiltonian_parameters, false);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_hamiltonian.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_hamiltonian;
    }

    // ALPHA separated evaluations
    SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_hamiltonian.outerSize(); ++i){
        for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_hamiltonian, i); it; ++it) {
            total_hamiltonian.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = oneElectronPartition(p, p, hamiltonian_parameters.get_g());
        const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, false);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_hamiltonian.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, hamiltonian_parameters.get_g());
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_hamiltonian.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    total_hamiltonian += this->calculateDiagonal(hamiltonian_parameters).asDiagonal();

    return total_hamiltonian;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the FCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
 *
 *  @return the action of the FCI Hamiltonian on the coefficient vector
 */
VectorX<double> FCI::matrixVectorProduct(const HamiltonianParameters<double>& hamiltonian_parameters, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("FCI::matrixVectorProduct(HamiltonianParameters<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    for (size_t p = 0; p<K; p++) {

        const auto& P = oneElectronPartition(p, p, hamiltonian_parameters.get_g());
        const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, false);

        // sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * xmap * this->alpha_couplings[p*(K+K+1-p)/2];
        for (size_t q = p + 1; q<K; q++) {

            const auto& P = oneElectronPartition(p, q, hamiltonian_parameters.get_g());
            const auto& beta_two_electron_intermediate = fock_space_beta.EvaluateOperatorSparse(P, true);

            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * xmap * this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
        }
    }

    auto beta_hamiltonian = fock_space_beta.EvaluateOperatorSparse(hamiltonian_parameters, false);
    auto alpha_hamiltonian = fock_space_alpha.EvaluateOperatorSparse(hamiltonian_parameters, false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian
 */
VectorX<double> FCI::calculateDiagonal(const HamiltonianParameters<double>& hamiltonian_parameters) const {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("FCI::calculateDiagonal(HamiltonianParameters<double>): Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto dim = fock_space.get_dimension();

    // Diagonal contributions
    VectorX<double> diagonal =  VectorX<double>::Zero(dim);

    auto k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();

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
                        diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += hamiltonian_parameters.get_g()(p, p, q, q);
                    }
                }  // q loop
            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
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



}  // namespace GQCP
