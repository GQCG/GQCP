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


#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"


namespace GQCP {
namespace CIEnvironment {


/**
 *  Create an environment suitable for solving dense CI eigenvalue problems for the given Hamiltonian and ONV basis.
 * 
 *  @tparam Hamiltonian             The type of Hamiltonian whose eigenproblem is trying to be solved.
 *  @tparam ONVBasis                The type of ONV basis in which the Hamiltonian should be represented.
 * 
 *  @param hamiltonian              A second-quantized Hamiltonian expressed in an orthonormal orbital basis.
 *  @param onv_basis                An ONV basis that spans a Fock (sub)space in which the Hamiltonian eigenproblem should be solved.
 * 
 *  @return An `EigenproblemEnvironment` initialized suitable for solving dense CI eigenvalue problems for the given Hamiltonian and ONV basis.
 */
template <typename Hamiltonian, typename ONVBasis>
EigenproblemEnvironment Dense(const Hamiltonian& hamiltonian, const ONVBasis& onv_basis) {

    // Determine the dense matrix representation of the Hamiltonian in the given ONV basis, and supply it to an `EigenproblemEnvironment`.
    const auto H = onv_basis.evaluateOperatorDense(hamiltonian);
    return EigenproblemEnvironment::Dense(H);
}


/**
 *  Create an environment suitable for solving iterative CI eigenvalue problems for the given Hamiltonian and ONV basis.
 * 
 *  @tparam Hamiltonian             The type of Hamiltonian whose eigenproblem is trying to be solved.
 *  @tparam ONVBasis                The type of ONV basis in which the Hamiltonian should be represented.
 * 
 *  @param hamiltonian              A second-quantized Hamiltonian expressed in an orthonormal orbital basis.
 *  @param onv_basis                An ONV basis that spans a Fock (sub)space in which the Hamiltonian eigenproblem should be solved.
 *  @param V                        A matrix of initial guess vectors, where each column of the matrix is an initial guess vector.
 * 
 *  @return An `EigenproblemEnvironment` initialized suitable for solving iterative CI eigenvalue problems for the given Hamiltonian and ONV basis.
 */
template <typename Hamiltonian, typename ONVBasis>
EigenproblemEnvironment Iterative(const Hamiltonian& hamiltonian, const ONVBasis& onv_basis, const MatrixX<double>& V) {

    // Determine the diagonal of the Hamiltonian matrix representation, and supply a matrix-vector product function to the `EigenproblemEnvironment`.
    const auto diagonal = onv_basis.evaluateOperatorDiagonal(hamiltonian);
    const auto matvec_function = [&hamiltonian, &onv_basis](const VectorX<double>& x) { return onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, x); };

    return EigenproblemEnvironment::Iterative(matvec_function, diagonal, V);
}


}  // namespace CIEnvironment
}  // namespace GQCP
