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
    size_t K;    // The number of spatial orbitals.
    size_t N_P;  // The number of electron pairs.


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The default constructor.
     */
    SeniorityZeroONVBasis() = default;

    /**
     *  @param K            The number of spatial orbitals.
     *  @param N_P          The number of electron pairs.
     */
    SeniorityZeroONVBasis(const size_t K, const size_t N_P);


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of electron pairs that this ONV basis is related to.
     */
    size_t numberOfElectronPairs() const { return this->N_P; }

    /**
     *  @return The number of spatial orbitals that this ONV basis is related to.
     */
    size_t numberOfSpatialOrbitals() const { return this->K; }

    /**
     *  Calculate the dimension of a seniority-zero ONV basis with a given number of spatial orbitals and number of electron pairs.
     * 
     *  @param K            The number of spatial orbitals.
     *  @param N_P          The number of electron pairs.
     * 
     *  @return The dimension of a seniority-zero ONV basis.
     */
    static size_t calculateDimension(const size_t K, const size_t N_P);

    /**
     *  @return The dimension of this ONV basis.
     */
    size_t dimension() const { return SeniorityZeroONVBasis::calculateDimension(this->numberOfSpatialOrbitals(), this->numberOfElectronPairs()); }


    /*
     *  MARK: Proxies
     */

    /**
     *  @return A spin-unresolved ONV basis that behaves analogously (with respect to a doubly-occupied situation) as this seniority-zero ONV basis.
     */
    SpinUnresolvedONVBasis proxy() const { return SpinUnresolvedONVBasis(this->K, this->N_P); }


    /*
     *  MARK: Iterations
     */

    /**
     *  Iterate over every (proxy) spin-(un)resolved ONV in this seniority-zero ONV basis and apply the given callback.
     * 
     *  @param callback             The function to be called on every step during the iteration over all the ONVs. The arguments of the callback are the ONV and its address in this ONV basis.
     */
    void forEach(const std::function<void(const SpinUnresolvedONV&, const size_t)>& callback) const;


    /*
     *  MARK: Dense restricted operator evaluations
     */

    /**
     *  Calculate the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return A dense matrix represention of the Hamiltonian.
     */
    SquareMatrix<double> evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Diagonal restricted operator evaluations
     */

    /**
     *  Calculate the diagonal of the matrix representation of a restricted one-electron operator in this ONV basis.
     *
     *  @param f_op             A restricted one-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the one-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const;

    /**
     *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
     *
     *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the two-electron operator.
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g_op) const;

    /**
     *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *
     *  @return The diagonal of the dense matrix represention of the Hamiltonian.
     */
    VectorX<double> evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const;


    /*
     *  MARK: Restricted matrix-vector product evaluations
     */

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted one-electron operator with the given coefficient vector.
     *
     *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const;

    /**
     *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
     *
     *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
     *  @param x                The coefficient vector of a linear expansion.
     *
     *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;
};


}  // namespace GQCP
