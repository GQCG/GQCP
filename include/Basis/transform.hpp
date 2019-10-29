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


#include "Basis/RSpinorBasis.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"



namespace GQCP {


/**
 *  Transform both the single-particle basis and the Hamiltonian to another basis using the given transformation matrix
 * 
 *  @tparam ShellType                   the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param sp_basis                     the single-particle basis
 *  @param sq_hamiltonian               the Hamiltonian
 *  @param T                            the transformation matrix
 */
template <typename TransformationScalar, typename ShellType>
void basisTransform(RSpinorBasis<TransformationScalar, ShellType>& sp_basis, SQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& T) {

    sp_basis.transform(T);
    sq_hamiltonian.transform(T);
}


/**
 *  Transform both the single-particle basis and the one-electron operator to another basis using the given transformation matrix
 * 
 *  @tparam ShellType                   the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param sp_basis                     the single-particle basis
 *  @param one_op                       the one-electron operator
 *  @param T                            the transformation matrix
 */
template <typename TransformationScalar, typename ShellType, size_t Components>
void basisTransform(RSpinorBasis<TransformationScalar, ShellType>& sp_basis, SQOneElectronOperator<TransformationScalar, Components>& one_op, const TransformationMatrix<TransformationScalar>& T) {

    sp_basis.transform(T);
    one_op.transform(T);
}


/**
 *  Rotate both the single-particle basis and the Hamiltonian to another basis using the given unitary transformation matrix
 * 
 *  @tparam ShellType                   the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param sp_basis                     the single-particle basis
 *  @param sq_hamiltonian               the Hamiltonian
 *  @param U                            the unitary transformation matrix
 */
template <typename TransformationScalar, typename ShellType>
void basisRotate(RSpinorBasis<TransformationScalar, ShellType>& sp_basis, SQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& U) {

    sp_basis.rotate(U);
    sq_hamiltonian.rotate(U);
}


/**
 *  Rotate both the single-particle basis and the Hamiltonian to another basis using the given Jacobi-rotation parameters
 * 
 *  @tparam ShellType                       the type of shell that the scalar basis contains
 * 
 *  @param sp_basis                         the single-particle basis
 *  @param sq_hamiltonian                   the Hamiltonian
 *  @param jacobi_rotation_parameters       the Jacobi-rotation parameters
 */
template <typename ShellType>
void basisRotate(RSpinorBasis<double, ShellType>& sp_basis, SQHamiltonian<double>& sq_hamiltonian, const JacobiRotationParameters& jacobi_rotation_parameters) {

    sp_basis.rotate(jacobi_rotation_parameters);
    sq_hamiltonian.rotate(jacobi_rotation_parameters);
}


}  // namespace GQCP
