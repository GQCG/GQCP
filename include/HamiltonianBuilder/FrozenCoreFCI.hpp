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


#include "FockSpace/FrozenProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/FrozenCoreCI.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A class capable of generating the matrix representation of the frozen core FCI Hamiltonian
 */
class FrozenCoreFCI : public FrozenCoreCI {
private:
    FrozenProductFockSpace fock_space;  // contains both the frozen alpha and beta Fock space

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the frozen product Fock space
     */
    explicit FrozenCoreFCI(const FrozenProductFockSpace& fock_space);


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // PUBLIC METHODS UNRESTRICTED METHODS
    using FrozenCoreCI::constructHamiltonian;
    using FrozenCoreCI::matrixVectorProduct;
    using FrozenCoreCI::calculateDiagonal;

    /**
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the FCI Hamiltonian matrix
     */
    SquareMatrix<double> constructHamiltonian(const USQHamiltonian<double>& usq_hamiltonian) const;

    /**
     *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param x                            the vector upon which the FCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
     *
     *  @return the action of the FCI Hamiltonian on the coefficient vector
     */
    VectorX<double> matrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    VectorX<double> calculateDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const;
};


}  // namespace GQCP
