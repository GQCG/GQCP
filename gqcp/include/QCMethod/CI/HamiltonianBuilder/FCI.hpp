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


#include "ONVBasis/ProductONVBasis.hpp"
#include "QCMethod/CI/HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {

/**
 *  A HamiltonianBuilder for FCI: it builds the matrix representation of the FCI Hamiltonian in the full alpha and beta product ONV basis
 */
class FCI : public HamiltonianBuilder {
private:
    ProductONVBasis fock_space;  // fock space containing the alpha and beta ONV basis


public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product ONV basis
     */
    explicit FCI(const ProductONVBasis& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN GETTERS
    const BaseONVBasis* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the FCI Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the FCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
     *
     *  @return the action of the FCI Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    VectorX<double> calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;
};


}  // namespace GQCP
