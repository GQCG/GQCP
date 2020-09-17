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


#include "Mathematical/Representation/SquareMatrix.hpp"
#include "ONVBasis/ONVBasisType.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <Eigen/Sparse>

#include <iostream>
#include <memory>


namespace GQCP {


/**
 *  A base class for the representation of an ONV basis
 */
class BaseONVBasis {
protected:
    size_t M;    // number of orbitals
    size_t dim;  // dimension of the ONV basis


public:
    // CONSTRUCTORS

    /**
     *  @param M        the number of orbitals
     *  @param dim      the dimension of the ONV basis
     */
    BaseONVBasis(const size_t M, const size_t dim);

    /**
     *  The default constructor.
     */
    BaseONVBasis() = default;


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    virtual ~BaseONVBasis() = default;


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @return the type of this ONV basis
     */
    virtual ONVBasisType type() const = 0;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const = 0;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const = 0;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param sq_hamiltonian              the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const = 0;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, const bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, const bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values              bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const = 0;


    // PUBLIC METHODS

    /**
     *  @return the dimension of this spin-unresolved ONV basis
     */
    size_t dimension() const { return this->dim; }

    /**
     *  @return the number of orbitals that this ONV basis describes
     */
    size_t numberOfOrbitals() const { return this->M; }
};


}  // namespace GQCP
