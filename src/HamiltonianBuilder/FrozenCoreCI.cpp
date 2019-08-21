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
#include "HamiltonianBuilder/FrozenCoreCI.hpp"

#include "FockSpace/BaseFrozenCoreFockSpace.hpp"
#include "Utilities/linalg.hpp"

#include <utility>


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param hamiltonian_builder           shared pointer to active (non-frozen core) Hamiltonian builder
 *  @param X                             the number of frozen orbitals
 */
FrozenCoreCI::FrozenCoreCI(std::shared_ptr<GQCP::HamiltonianBuilder> hamiltonian_builder, size_t X) :
    HamiltonianBuilder(),
    active_hamiltonian_builder (std::move(hamiltonian_builder)),
    X (X)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the frozen core Hamiltonian matrix
 */
SquareMatrix<double> FrozenCoreCI::constructHamiltonian(const HamiltonianParameters<double>& ham_par) const {

    // Freeze Hamiltonian parameters
    HamiltonianParameters<double> frozen_ham_par =  BaseFrozenCoreFockSpace::freezeOperator(ham_par, X);

    // calculate Hamiltonian matrix through conventional CI
    SquareMatrix<double> total_hamiltonian = this->active_hamiltonian_builder->constructHamiltonian(frozen_ham_par);

    // diagonal correction
    VectorX<double> diagonal = VectorX<double>::Ones(this->get_fock_space()->get_dimension());
    auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(ham_par, this->X, active_hamiltonian_builder->get_fock_space()->get_dimension());
    total_hamiltonian += frozen_core_diagonal.asDiagonal();

    return total_hamiltonian;
}


/**
 *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x            the vector upon which the Hamiltonian acts
 *  @param diagonal     the diagonal of the Hamiltonian matrix
 *
 *  @return the action of the frozen core Hamiltonian on the coefficient vector
 */
VectorX<double> FrozenCoreCI::matrixVectorProduct(const HamiltonianParameters<double>& ham_par, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    HamiltonianParameters<double> frozen_ham_par =  BaseFrozenCoreFockSpace::freezeOperator(ham_par, X);

    // perform matvec in the active space with "frozen" Hamiltonian parameters
    return this->active_hamiltonian_builder->matrixVectorProduct(frozen_ham_par, x, diagonal);
}


/**
 *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
 */
VectorX<double> FrozenCoreCI::calculateDiagonal(const HamiltonianParameters<double>& ham_par) const {

    HamiltonianParameters<double> frozen_ham_par =  BaseFrozenCoreFockSpace::freezeOperator(ham_par, this->X);

    // calculate diagonal in the active space with the "frozen" Hamiltonian parameters
    VectorX<double> diagonal = this->active_hamiltonian_builder->calculateDiagonal(frozen_ham_par);

    // calculate diagonal for the frozen orbitals
    auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(ham_par, this->X, active_hamiltonian_builder->get_fock_space()->get_dimension());

    return diagonal + frozen_core_diagonal;
}


}  // namespace GQCP

