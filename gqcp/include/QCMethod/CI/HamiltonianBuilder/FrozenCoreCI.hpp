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
#pragma once


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {


/**
 *  (base) Class implementing general functions related to frozen core CI
 */
class FrozenCoreCI: public HamiltonianBuilder {
protected:
    size_t X;                                                        // number of frozen orbitals/electrons
    std::shared_ptr<HamiltonianBuilder> active_hamiltonian_builder;  // non-frozen core Hamiltonian builder performing the HamiltonianBuilder interface in the active space with the frozen Hamiltonian

public:
    // CONSTRUCTORS
    /**
     *  @param hamiltonian_builder           shared pointer to active (non-frozen core) Hamiltonian builder
     *  @param X                             the number of frozen orbitals
     */
    FrozenCoreCI(std::shared_ptr<HamiltonianBuilder> hamiltonian_builder, size_t X);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the frozen core Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param x                    the vector upon which the Hamiltonian acts
     *  @param diagonal             the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the frozen core Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
     */
    VectorX<double> calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;
};


}  // namespace GQCP
