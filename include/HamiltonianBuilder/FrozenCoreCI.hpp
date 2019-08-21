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


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/**
 *  (base) Class implementing general functions related to frozen core CI
 */
class FrozenCoreCI : public HamiltonianBuilder {
protected:
    size_t X;  // number of frozen orbitals/electrons
    std::shared_ptr<HamiltonianBuilder> active_hamiltonian_builder;  // non-frozen core Hamiltonian builder performing the HamiltonianBuilder interface in the active space with the frozen Hamiltonian parameters

public:
    // CONSTRUCTORS
    /**
     *  @param hamiltonian_builder           shared pointer to active (non-frozen core) Hamiltonian builder
     *  @param X                             the number of frozen orbitals
     */
    FrozenCoreCI(std::shared_ptr<HamiltonianBuilder> hamiltonian_builder, size_t X);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the frozen core Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const HamiltonianParameters<double>& ham_par) const override;

    /**
     *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x            the vector upon which the Hamiltonian acts
     *  @param diagonal     the diagonal of the Hamiltonian matrix
     *
     *  @return the action of the frozen core Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const HamiltonianParameters<double>& ham_par, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

    /**
     *  @param ham_par      the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
     */
    VectorX<double> calculateDiagonal(const HamiltonianParameters<double>& ham_par) const override;
};


}  // namespace GQCP
