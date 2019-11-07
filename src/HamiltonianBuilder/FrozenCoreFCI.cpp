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
#include "HamiltonianBuilder/FrozenCoreFCI.hpp"

#include "HamiltonianBuilder/FCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the frozen product Fock space
 */
FrozenCoreFCI::FrozenCoreFCI(const FrozenProductFockSpace& fock_space) :
    FrozenCoreCI(std::make_shared<FCI>(fock_space.get_active_product_fock_space()), fock_space.get_number_of_frozen_orbitals()),
    fock_space (fock_space)
{}

/*
 *  UNRESTRICTED METHODS
 */

/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an unrestricted orthonormal basis
 *
 *  @return the frozen core Hamiltonian matrix
 */
SquareMatrix<double> FrozenCoreFCI::constructHamiltonian(const USQHamiltonian<double>& sq_hamiltonian) const {

    // Freeze the Hamiltonian
    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, X);
    const auto& fci = static_cast<const FCI&>(*this->active_hamiltonian_builder);

    // calculate Hamiltonian matrix through conventional CI
    SquareMatrix<double> total_hamiltonian = fci.constructHamiltonian(frozen_ham_par);

    // diagonal correction
    VectorX<double> diagonal = VectorX<double>::Ones(this->fock_space.get_dimension());
    auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian, this->X, fci.get_fock_space()->get_dimension());
    total_hamiltonian += frozen_core_diagonal.asDiagonal();

    return total_hamiltonian;
}


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an unrestricted orthonormal basis
 *  @param x                    the vector upon which the Hamiltonian acts
 *  @param diagonal             the diagonal of the Hamiltonian matrix
 *
 *  @return the action of the frozen core Hamiltonian on the coefficient vector
 */
VectorX<double> FrozenCoreFCI::matrixVectorProduct(const USQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, X);
    const auto& fci = static_cast<const FCI&>(*this->active_hamiltonian_builder);


    // Perform the matvec in the active space with the "frozen"
    return fci.matrixVectorProduct(frozen_ham_par, x, diagonal);
}


/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an unrestricted orthonormal basis
 *
 *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
 */
VectorX<double> FrozenCoreFCI::calculateDiagonal(const USQHamiltonian<double>& sq_hamiltonian) const {

    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, this->X);
    const auto& fci = static_cast<const FCI&>(*this->active_hamiltonian_builder);

    // Calculate the diagonal in the active space with the "frozen" Hamiltonian
    VectorX<double> diagonal = fci.calculateDiagonal(frozen_ham_par);

    // Calculate the diagonal for the frozen orbitals
    auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian, this->X, fci.get_fock_space()->get_dimension());

    return diagonal + frozen_core_diagonal;
}


}  // namespace GQCP
