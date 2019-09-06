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
{}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the FCI Hamiltonian matrix
 */
SquareMatrix<double> FCI::constructHamiltonian(const SQHamiltonian<double>& hamiltonian_parameters) const {

    return this->fock_space.evaluateOperatorDense(hamiltonian_parameters, true);
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the FCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
 *
 *  @return the action of the FCI Hamiltonian on the coefficient vector
 */
VectorX<double> FCI::matrixVectorProduct(const SQHamiltonian<double>& hamiltonian_parameters, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = hamiltonian_parameters.core().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("FCI::matrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    const auto& alpha_couplings = this->fock_space.get_alpha_couplings();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    for (size_t p = 0; p<K; p++) {

        const auto& P = this->fock_space.oneElectronPartition(p, p, hamiltonian_parameters.twoElectron());
        const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, false);

        // sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2]);
        for (size_t q = p + 1; q<K; q++) {

            const auto& P = this->fock_space.oneElectronPartition(p, q, hamiltonian_parameters.twoElectron());
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2 + q - p]);
        }
    }

    auto beta_hamiltonian = fock_space_beta.evaluateOperatorSparse(hamiltonian_parameters, false);
    auto alpha_hamiltonian = fock_space_alpha.evaluateOperatorSparse(hamiltonian_parameters, false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian
 */
VectorX<double> FCI::calculateDiagonal(const SQHamiltonian<double>& hamiltonian_parameters) const {

   return this->fock_space.evaluateOperatorDiagonal(hamiltonian_parameters);
}



}  // namespace GQCP
