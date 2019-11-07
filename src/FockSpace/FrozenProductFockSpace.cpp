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
#include "FockSpace/FrozenProductFockSpace.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the total number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the total number of alpha electrons
 *  @param N_beta       the total number of beta electrons
 *  @param X            the number of frozen orbitals and electrons (equal for alpha and beta)
 */
FrozenProductFockSpace::FrozenProductFockSpace(size_t K, size_t N_alpha, size_t N_beta, size_t X) :
        BaseFrozenCoreFockSpace(std::make_shared<ProductFockSpace>(ProductFockSpace(K-X, N_alpha-X, N_beta-X)), X),
        frozen_fock_space_alpha (FrozenFockSpace(K, N_alpha, X)),
        frozen_fock_space_beta (FrozenFockSpace(K, N_beta, X)),
        active_product_fock_space (ProductFockSpace(K-X, N_alpha-X, N_beta-X)),
        X (X)
{}


/**
 *  @param fock_space       (to be frozen) product Fock space
 *  @param X                the number of frozen orbitals and electrons (equal for alpha and beta)
 */
FrozenProductFockSpace::FrozenProductFockSpace(const ProductFockSpace& fock_space, size_t X) :
    FrozenProductFockSpace(fock_space.get_K(), fock_space.get_N_alpha(), fock_space.get_N_beta(), X)
{}


/*
 *  UNRESTRICTED
 */

/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param diagonal_values                bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> FrozenProductFockSpace::evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const {
    // Freeze the operators
    const auto frozen_usq_hamiltonian = BaseFrozenCoreFockSpace::freezeOperator(usq_hamiltonian, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_product_fock_space.evaluateOperatorDense(frozen_usq_hamiltonian, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(usq_hamiltonian, this->X, this->active_product_fock_space.get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}

/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param usq_hamiltonian_alpha          the Hamiltonian expressed in an unrestricted orthonormal basis 
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> FrozenProductFockSpace::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {
    const auto frozen_usq_hamiltonian = BaseFrozenCoreFockSpace::freezeOperator(usq_hamiltonian, this->X);

    // Calculate diagonal in the active space with the "frozen" Hamiltonian
    const auto diagonal = this->active_product_fock_space.evaluateOperatorDiagonal(frozen_usq_hamiltonian);

    // Calculate diagonal for the frozen orbitals
    const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(usq_hamiltonian, this->X, this->active_product_fock_space.get_dimension());

    return diagonal + frozen_core_diagonal;

}


}  // namespace GQCP
