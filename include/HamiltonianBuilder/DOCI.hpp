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


#include "FockSpace/FockSpace.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"

#include <memory>


namespace GQCP {


/**
 *  A HamiltonianBuilder for DOCI: it builds the matrix representation of the DOCI Hamiltonian, in a Fock space where orbitals are either doubly occupied or unoccupied.
 */
class DOCI : public HamiltonianBuilder {
private:
    FockSpace fock_space;  // both the alpha and beta Fock space


public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the full Fock space, identical for alpha and beta
     */
    explicit DOCI(const FockSpace& fock_space);


    // DESTRUCTOR
    ~DOCI() = default;


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the DOCI Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the DOCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the DOCI Hamiltonian matrix
     *
     *  @return the action of the DOCI Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the diagonal of the matrix representation of the DOCI Hamiltonian
     */
    VectorX<double> calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;
};


}  // namespace GQCP
