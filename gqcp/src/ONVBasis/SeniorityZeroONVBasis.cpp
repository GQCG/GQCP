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
#include "ONVBasis/SeniorityZeroONVBasis.hpp"

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 */
SeniorityZeroONVBasis::SeniorityZeroONVBasis(const size_t K, const size_t N_P) :
    K (K),
    N_P (N_P),
    dim (SeniorityZeroONVBasis::calculateDimension(K, N_P))
{}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 * 
 *  @return the dimension of a seniority-zero ONV basis with the given number of spatial orbitals and electron pairs
 */
size_t SeniorityZeroONVBasis::calculateDimension(const size_t K, const size_t N_P) {

    return SpinUnresolvedONVBasis::calculateDimension(K, N_P);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Evaluate the diagonal of the operator
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
 */
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const SpinResolvedSelectedONVBasis selected_onv_basis {*this};
    return selected_onv_basis.evaluateOperatorDiagonal(one_op);
}


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the ONV basis
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of the ONV basis
 */
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

    const SpinResolvedSelectedONVBasis selected_onv_basis {*this};
    return selected_onv_basis.evaluateOperatorMatrixVectorProduct(one_op, x, diagonal);
}


/**
 *  @return a coefficient vector that describes the expansion coefficients of the Hartree-Fock wave function (i.e. the single Slater determinant with the lowest energy)
 */
VectorX<double> SeniorityZeroONVBasis::hartreeFockExpansion() const {

    VectorX<double> coefficients = VectorX<double>::Zero(this->dim);
    coefficients(0) = 1;  // the lowest-energy SSD has position 0: this is conventially determined by the ONV basis
    return coefficients;
}


/**
 *  @return a coefficient vector that describes the expansion coefficients of a random, normalized linear expansion
 */
VectorX<double> SeniorityZeroONVBasis::randomExpansion() const {

    VectorX<double> coefficients = VectorX<double>::Random(this->dim);
    coefficients.normalize();
    return coefficients;
}


}  // namespace GQCP
