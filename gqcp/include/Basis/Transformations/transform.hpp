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


#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/TransformationMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QuantumChemical/Spin.hpp"


namespace GQCP {


/**
 *  Transform both the spinor basis and the Hamiltonian to another basis using the given transformation matrix
 * 
 *  @tparam Shell                       the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param sq_hamiltonian               the Hamiltonian
 *  @param T                            the transformation matrix
 */
template <typename TransformationScalar, typename Shell>
void basisTransform(RSpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, RSQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& T) {

    spinor_basis.transform(T);
    sq_hamiltonian.transform(T);
}

template <typename TransformationScalar, typename Shell>
void basisTransform(RSpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, RSQHamiltonian<TransformationScalar>& sq_hamiltonian, const RTransformation<TransformationScalar>& T) {

    spinor_basis.transform(T);
    sq_hamiltonian.transform(T);
}

/**
 *  Transform both the spinor basis and the one-electron operator to another basis using the given transformation matrix
 * 
 *  @tparam Shell                       the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param one_op                       the one-electron operator
 *  @param T                            the transformation matrix
 */
template <typename TransformationScalar, typename Shell, typename Vectorizer>
void basisTransform(RSpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, RSQOneElectronOperator<TransformationScalar, Vectorizer>& one_op, const TransformationMatrix<TransformationScalar>& T) {

    spinor_basis.transform(T);
    one_op.transform(T);
}


/**
 *  Rotate both the spinor basis and the Hamiltonian to another basis using the given unitary transformation matrix
 * 
 *  @tparam Shell                       the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param sq_hamiltonian               the Hamiltonian
 *  @param U                            the unitary transformation matrix
 */
template <typename TransformationScalar, typename Shell>
void basisRotate(RSpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, RSQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& U) {

    spinor_basis.rotate(U);
    sq_hamiltonian.rotate(U);
}

template <typename TransformationScalar, typename Shell>
void basisRotate(RSpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, RSQHamiltonian<TransformationScalar>& sq_hamiltonian, const RTransformation<TransformationScalar>& U) {

    spinor_basis.rotate(U);
    sq_hamiltonian.rotate(U);
}


/**
 *  Rotate both the spinor basis and the Hamiltonian to another basis using the given Jacobi rotation.
 * 
 *  @tparam Shell                           the type of shell that the scalar basis contains
 * 
 *  @param spinor_basis                     the spinor basis
 *  @param sq_hamiltonian                   the Hamiltonian
 *  @param jacobi_rotation                  The Jacobi rotation.
 */
template <typename Shell>
void basisRotate(RSpinOrbitalBasis<double, Shell>& spinor_basis, RSQHamiltonian<double>& sq_hamiltonian, const JacobiRotation& jacobi_rotation) {

    spinor_basis.rotate(jacobi_rotation);
    sq_hamiltonian.rotate(jacobi_rotation);
}


/**
 *  UNRESTRICTED
 */

/**
 *  Transform a single component (alpha or beta) of both the spinor basis and the Hamiltonian to another basis using the given transformation matrix
 * 
 *  @tparam ShellType                   the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param usq_hamiltonian              the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param T                            the transformation matrix for one of the spin components
 *  @param component                    the spin component
 */
// template <typename TransformationScalar, typename ShellType>
// void basisTransform(USpinOrbitalBasis<TransformationScalar, ShellType>& spinor_basis, USQHamiltonian<TransformationScalar>& usq_hamiltonian, const TransformationMatrix<TransformationScalar>& T, const Spin& component) {
//     spinor_basis.transform(T, component);
//     usq_hamiltonian.transform(T, component);
// }


/**
 *  Transform both the spinor basis and the Hamiltonian to another basis using the given transformation matrix
 * 
 *  @tparam ShellType                   the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param usq_hamiltonian              the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param T                            the transformation matrix for both of the spin components
 */
// template <typename TransformationScalar, typename ShellType>
// void basisTransform(USpinOrbitalBasis<TransformationScalar, ShellType>& spinor_basis, USQHamiltonian<TransformationScalar>& usq_hamiltonian, const TransformationMatrix<TransformationScalar>& T) {
//     spinor_basis.transform(T);
//     usq_hamiltonian.transform(T);
// }
// template <typename TransformationScalar, typename ShellType>
// void basisTransform(USpinOrbitalBasis<TransformationScalar, ShellType>& spinor_basis, USQHamiltonian<TransformationScalar>& usq_hamiltonian, const UTransformationComponent<TransformationScalar>& T) {
//     spinor_basis.transform(T);
//     usq_hamiltonian.transform(T);
// }

/**
 *  Rotate a single component (alpha or beta) of both the spinor basis and the Hamiltonian to another basis using the given unitary transformation matrix
 *  
 *  @tparam Shell                       the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param usq_hamiltonian              the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param U                            the unitary transformation matrix for the given spin component
 *  @param component                    the spin component
 */
// template <typename TransformationScalar, typename Shell>
// void basisRotate(USpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, USQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& U, const Spin& component) {

//     spinor_basis.rotate(U, component);
//     sq_hamiltonian.rotate(U, component);
// }


/**
 *  Rotate both components of the both the spinor basis and the Hamiltonian to another basis using the given unitary transformation matrix
 *  
 *  @tparam Shell                       the type of shell that the scalar basis contains
 *  @tparam TransformationScalar        the scalar type of the transformation matrix
 * 
 *  @param spinor_basis                 the spinor basis
 *  @param usq_hamiltonian              the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param U                            the unitary transformation matrix for both of the spin components
 */
// template <typename TransformationScalar, typename Shell>
// void basisRotate(USpinOrbitalBasis<TransformationScalar, Shell>& spinor_basis, USQHamiltonian<TransformationScalar>& sq_hamiltonian, const TransformationMatrix<TransformationScalar>& U) {

//     spinor_basis.rotate(U);
//     sq_hamiltonian.rotate(U);
// }


/**
 *  Rotate one spin component of both the spinor basis and the Hamiltonian to another basis using the given Jacobi-rotation parameters
 *  
 *  @tparam Shell                           the type of shell that the scalar basis contains
 * 
 *  @param spinor_basis                     the spinor basis
 *  @param usq_hamiltonian                  the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param jacobi_rotation_parameters       the Jacobi-rotation parameters for one of the spin components
 *  @param component                        the spin component
 */
// template <typename Shell>
// void basisRotate(USpinOrbitalBasis<double, Shell>& spinor_basis, USQHamiltonian<double>& sq_hamiltonian, const JacobiRotation& jacobi_rotation_parameters, const Spin& component) {

//     spinor_basis.rotate(jacobi_rotation_parameters, component);
//     sq_hamiltonian.rotate(jacobi_rotation_parameters, component);
// }


/**
 *  Rotate both components of both the spinor basis and the Hamiltonian to another basis using the given Jacobi-rotation parameters
 * 
 *  @tparam Shell                           the type of shell that the scalar basis contains
 * 
 *  @param spinor_basis                     the single-particle basis for the beta component
 *  @param usq_hamiltonian                  the Hamiltonian, expressed in an unrestricted spinor basis
 *  @param jacobi_rotation_parameters       the Jacobi-rotation parameters
 */
// template <typename Shell>
// void basisRotate(USpinOrbitalBasis<double, Shell>& spinor_basis, USQHamiltonian<double>& sq_hamiltonian, const JacobiRotation& jacobi_rotation_parameters) {

//     spinor_basis.rotate(jacobi_rotation_parameters);
//     sq_hamiltonian.rotate(jacobi_rotation_parameters);
// }


}  // namespace GQCP
