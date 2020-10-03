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

#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreCI.hpp"

#include "ONVBasis/BaseFrozenCoreONVBasis.hpp"
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
FrozenCoreCI::FrozenCoreCI(const std::shared_ptr<GQCP::HamiltonianBuilder> hamiltonian_builder, const size_t X) :
    HamiltonianBuilder(),
    active_hamiltonian_builder {std::move(hamiltonian_builder)},
    X {X} {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
 */
VectorX<double> FrozenCoreCI::calculateDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const {

    RSQHamiltonian<double> frozen_ham_par = BaseFrozenCoreONVBasis::freezeOperator(sq_hamiltonian, this->X);

    // Calculate the diagonal in the active space with the "frozen" Hamiltonian
    VectorX<double> diagonal = this->active_hamiltonian_builder->calculateDiagonal(frozen_ham_par);

    // Calculate the diagonal for the frozen orbitals
    auto frozen_core_diagonal = BaseFrozenCoreONVBasis::frozenCoreDiagonal(sq_hamiltonian, this->X, active_hamiltonian_builder->onvBasis()->dimension());

    return diagonal + frozen_core_diagonal;
}


/**
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the frozen core Hamiltonian matrix
 */
SquareMatrix<double> FrozenCoreCI::constructHamiltonian(const RSQHamiltonian<double>& sq_hamiltonian) const {

    // Freeze the Hamiltonian
    RSQHamiltonian<double> frozen_ham_par = BaseFrozenCoreONVBasis::freezeOperator(sq_hamiltonian, X);

    // calculate Hamiltonian matrix through conventional CI
    SquareMatrix<double> total_hamiltonian = this->active_hamiltonian_builder->constructHamiltonian(frozen_ham_par);

    // diagonal correction
    VectorX<double> diagonal = VectorX<double>::Ones(this->onvBasis()->dimension());
    auto frozen_core_diagonal = BaseFrozenCoreONVBasis::frozenCoreDiagonal(sq_hamiltonian, this->X, active_hamiltonian_builder->onvBasis()->dimension());
    total_hamiltonian += frozen_core_diagonal.asDiagonal();

    return total_hamiltonian;
}


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param x                    the vector upon which the Hamiltonian acts
 *  @param diagonal             the diagonal of the Hamiltonian matrix
 *
 *  @return the action of the frozen core Hamiltonian on the coefficient vector
 */
VectorX<double> FrozenCoreCI::matrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    RSQHamiltonian<double> frozen_ham_par = BaseFrozenCoreONVBasis::freezeOperator(sq_hamiltonian, X);

    // Perform the matvec in the active space with the "frozen"
    return this->active_hamiltonian_builder->matrixVectorProduct(frozen_ham_par, x, diagonal);
}


}  // namespace GQCP
