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


#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "ONVBasis/SpinResolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedFrozenONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"


namespace GQCP {
namespace CIEnvironment {


/**
 *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
 *  @param onv_basis                        the full, spin-resolved ONV basis
 * 
 *  @return an environment suitable for solving Hubbard-related eigenvalue problems
 */
template <typename Scalar>
EigenproblemEnvironment Dense(const HubbardHamiltonian<Scalar>& hubbard_hamiltonian, const SpinResolvedONVBasis& onv_basis) {

    const Hubbard hubbard_builder (onv_basis);  // the 'HamiltonianBuilder'
    const auto H = hubbard_builder.constructHamiltonian(hubbard_hamiltonian);
    return EigenproblemEnvironment::Dense(H);
}


/**
 *  @param sq_hamiltonian               the general, second-quantized representation of the Hamiltonian
 *  @param onv_basis                    the full, spin-resolved ONV basis
 * 
 *  @return an environment suitable for solving spin-resolved FCI eigenvalue problems
 */
template <typename Scalar>
EigenproblemEnvironment Dense(const SQHamiltonian<Scalar>& sq_hamiltonian, const SpinResolvedONVBasis& onv_basis) {

    const FCI fci_builder (onv_basis);  // the 'HamiltonianBuilder'
    const auto H = fci_builder.constructHamiltonian(sq_hamiltonian);
    return EigenproblemEnvironment::Dense(H);
}


/**
 *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
 *  @param onv_basis                        the full, spin-resolved ONV basis
 *  @param V                                a matrix of initial guess vectors (each column of the matrix is an initial guess vector)
 * 
 *  @return an environment suitable for solving Hubbard-related eigenvalue problems
 */
template <typename Scalar>
EigenproblemEnvironment Iterative(const HubbardHamiltonian<Scalar>& hubbard_hamiltonian, const SpinResolvedONVBasis& onv_basis, const MatrixX<double>& V) {

    const Hubbard hubbard_builder (onv_basis);  // the 'HamiltonianBuilder'

    const auto diagonal = hubbard_builder.calculateDiagonal(hubbard_hamiltonian);
    const auto matvec_function = [diagonal, hubbard_builder, hubbard_hamiltonian] (const VectorX<double>& x) { return hubbard_builder.matrixVectorProduct(hubbard_hamiltonian, x, diagonal); };

    return EigenproblemEnvironment::Iterative(matvec_function, diagonal, V);
}


/**
 *  @param sq_hamiltonian               the general, second-quantized representation of the Hamiltonian
 *  @param onv_basis                    the full, spin-resolved ONV basis
 *  @param V                            a matrix of initial guess vectors (each column of the matrix is an initial guess vector)
 * 
 *  @return an environment suitable for solving spin-resolved FCI eigenvalue problems
 */
template <typename Scalar>
EigenproblemEnvironment Iterative(const SQHamiltonian<Scalar>& sq_hamiltonian, const SpinResolvedONVBasis& onv_basis, const MatrixX<double>& V) {

    const FCI fci_builder (onv_basis);  // the 'HamiltonianBuilder'

    const auto diagonal = fci_builder.calculateDiagonal(sq_hamiltonian);
    const auto matvec_function = [diagonal, fci_builder, sq_hamiltonian] (const VectorX<double>& x) { return fci_builder.matrixVectorProduct(sq_hamiltonian, x, diagonal); };

    return EigenproblemEnvironment::Iterative(matvec_function, diagonal, V);
}


}  // namespace CIEnvironment
}  // namespace GQCP
