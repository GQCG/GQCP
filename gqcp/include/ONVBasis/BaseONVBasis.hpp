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
    size_t K;    // number of spatial orbitals
    size_t dim;  // dimension of the ONV basis


    // PROTECTED CONSTRUCTORS
    BaseONVBasis() = default;
    /**
     *  @param K        the number of orbitals
     *  @param dim      the dimension of the ONV basis
     */
    BaseONVBasis(size_t K, size_t dim);


public:
    // NAMED CONSTRUCTORS

    /**
     *  Clones a derived BaseONVBasis instance to the heap memory
     *
     *  @param fock_space     reference to a derived BaseONVBasis instance to be cloned.
     *
     *  @return a shared pointer owning the heap-cloned ONV basis
     */
    static std::shared_ptr<BaseONVBasis> CloneToHeap(const BaseONVBasis& fock_space);


    // DESTRUCTOR
    virtual ~BaseONVBasis() = default;


    // GETTERS
    size_t get_dimension() const { return dim; }
    size_t get_K() const { return K; }
    virtual ONVBasisType get_type() const = 0;


    // PUBLIC METHODS

    /**
     *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
     */
    VectorX<double> hartreeFockExpansion() const;

    /**
     *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
     */
    VectorX<double> randomExpansion() const;

    /**
     *  @return a constant normalized coefficients vector (i.e. all the coefficients are equal)
     */
    VectorX<double> constantExpansion() const;


    // PUBLIC VIRTUAL METHODS

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values              bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const = 0;

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
};


}  // namespace GQCP
