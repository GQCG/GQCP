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

#pragma once


#include "ONVBasis/BaseONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <memory>
#include <utility>


namespace GQCP {


/**
 *  A base class whose derived classes are able to construct matrix representations of the Hamiltonian in an ONV basis
 *
 *  Derived classes should implement:
 *      - constructHamiltonian() which constructs the full Hamiltonian matrix in the given ONV basis
 *      - matrixVectorProduct() which gives the result of the action of the Hamiltonian on a given coefficient vector
 *      - calculateDiagonal() which gives the diagonal of the Hamiltonian matrix
 */
class HamiltonianBuilder {
public:
    // DESTRUCTOR
    virtual ~HamiltonianBuilder() = default;


    // PURE VIRTUAL GETTERS
    virtual const BaseONVBasis* get_fock_space() const = 0;


    // PURE VIRTUAL PUBLIC METHODS

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian matrix
     */
    virtual SquareMatrix<double> constructHamiltonian(const SQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the Hamiltonian on the coefficient vector
     */
    virtual VectorX<double> matrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const = 0;

    /**
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    virtual VectorX<double> calculateDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const = 0;
};


}  // namespace GQCP
