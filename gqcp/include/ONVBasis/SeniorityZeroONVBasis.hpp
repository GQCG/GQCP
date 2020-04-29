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


#include "Mathematical/Representation/Matrix.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <functional>


namespace GQCP {


/**
 *  A spin-resolved ONV basis that contains all seniority-zero (i.e. doubly-occupied) (spin-resolved) ONVs.
 */
class SeniorityZeroONVBasis {
private:
    size_t K;    // the number of spatial orbitals
    size_t N_P;  // the number of electron pairs

    size_t dim;  // the dimension of this ONV basis


public:
    // CONSTRUCTORS

    /**
     *  The default constructor.
     */
    SeniorityZeroONVBasis() = default;


    /**
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     */
    SeniorityZeroONVBasis(const size_t K, const size_t N_P);


    // STATIC PUBLIC METHODS

    /**
     *  @param K            the number of spatial orbitals
     *  @param N_P          the number of electron pairs
     * 
     *  @return the dimension of a seniority-zero ONV basis with the given number of spatial orbitals and electron pairs
     */
    static size_t calculateDimension(const size_t K, const size_t N_P);


    // PUBLIC METHODS

    /**
     *  @return the dimension of this ONV basis
     */
    size_t dimension() const { return this->dim; }

    /**
     *  Evaluate the diagonal of the matrix representation of a one-electron operator inside this seniority-zero ONV basis.
     *
     *  @param one_op               a one-electron operator expressed in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the one-electron operator in this seniority-zero ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const;

    /**
     *  Evaluate the diagonal of the matrix representation of a two-electron operator inside this seniority-zero ONV basis.
     *
     *  @param two_op               a two-electron operator expressed in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the two-electron operator in this seniority-zero ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const;

    /**
     *  Evaluate the diagonal of the matrix representation of a Hamiltonian inside this seniority-zero ONV basis.
     *
     *  @param sq_hamiltonian               a Hamiltonian expressed in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian in this seniority-zero ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const;

    /**
     *  Evaluate a one electron operator in a matrix vector product
     *
     *  @param one_op                       the one electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the ONV basis
     *
     *  @return the one electron operator's matrix vector product in a vector with the dimensions of the ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Iterate over every (proxy) spin-resolved ONV in this seniority-zero ONV basis and apply the given callback.
     * 
     *  @param callback             a function to be called on every step during the iteration over all the ONVs. The arguments of the callback are the ONV and its address in this ONV basis.
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const;

    /**
     *  @return a coefficient vector that describes the expansion coefficients of the Hartree-Fock wave function (i.e. the single Slater determinant with the lowest energy)
     */
    VectorX<double> hartreeFockExpansion() const;

    /**
     *  @return the number of electron pairs that this ONV basis is related to
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return the number of spatial orbitals that this ONV basis is related to
     */
    size_t numberOfSpatialOrbitals() const { return this->K; }

    /**
     *  @return a spin-unresolved ONV basis that behaves analogously (with respect to a doubly-occupied situation) as this seniority-zero ONV basis
     */
    SpinUnresolvedONVBasis proxy() const { return SpinUnresolvedONVBasis(this->K, this->N_P); }

    /**
     *  @return a coefficient vector that describes the expansion coefficients of a random, normalized linear expansion
     */
    VectorX<double> randomExpansion() const;
};


}  // namespace GQCP
