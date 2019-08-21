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


#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {


/**
 *  Hubbard builds a a Hubbard Hamiltonian matrix in the FCI Fock space
 *
 *  Hubbard distinguishes itself from FCI by explicitly implementing simplified Hamiltonian parameters:
 *      - for the one electron operators only inter-site interactions are considered
 *      - for the two electron operators only on-site (doubly occupied in-place) interactions are considered
 */
class Hubbard : public HamiltonianBuilder {
private:
    ProductFockSpace fock_space;  // fock space containing the alpha and beta Fock space
    

public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit Hubbard(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~Hubbard() = default;


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the Hubbard Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const HamiltonianParameters<double>& hamiltonian_parameters) const override;

    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the Hubbard Hamiltonian acts
     *  @param diagonal                     the diagonal of the Hubbard Hamiltonian matrix
     *
     *  @return the action of the Hubbard Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const HamiltonianParameters<double>& hamiltonian_parameters, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

    /**
     *  @param hamiltonian_parameters       the Hubbard Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hubbard Hamiltonian
     */
    VectorX<double> calculateDiagonal(const HamiltonianParameters<double>& hamiltonian_parameters) const override;
};



}  // namespace GQCP
