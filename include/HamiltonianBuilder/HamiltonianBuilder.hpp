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
#ifndef GQCP_HAMILTONIANBUILDER_HPP
#define GQCP_HAMILTONIANBUILDER_HPP


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "FockSpace/BaseFockSpace.hpp"

#include <memory>
#include <utility>



namespace GQCP {


/**
 *  A base class whose derived classes are able to construct matrix representations of the Hamiltonian in a Fock space
 *
 *  Derived classes should implement:
 *      - constructHamiltonian() which constructs the full Hamiltonian matrix in the given Fock space
 *      - matrixVectorProduct() which gives the result of the action of the Hamiltonian on a given coefficient vector
 *      - calculateDiagonal() which gives the diagonal of the Hamiltonian matrix
 */
class HamiltonianBuilder {
public:
    // DESTRUCTOR
    virtual ~HamiltonianBuilder() = default;


    // PURE VIRTUAL GETTERS
    virtual const BaseFockSpace* get_fock_space() const = 0;


    // PURE VIRTUAL PUBLIC METHODS

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the Hamiltonian matrix
     */
    virtual SquareMatrix<double> constructHamiltonian(const HamiltonianParameters<double>& hamiltonian_parameters) const = 0;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the Hamiltonian on the coefficient vector
     */
    virtual VectorX<double> matrixVectorProduct(const HamiltonianParameters<double>& hamiltonian_parameters, const VectorX<double>& x, const VectorX<double>& diagonal) const = 0;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    virtual VectorX<double> calculateDiagonal(const HamiltonianParameters<double>& hamiltonian_parameters) const = 0;
};


}  // namespace GQCP



#endif  // GQCP_HAMILTONIANBUILDER_HPP
