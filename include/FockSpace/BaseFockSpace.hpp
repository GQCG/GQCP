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
#ifndef GQCP_BASEFOCKSPACE_HPP
#define GQCP_BASEFOCKSPACE_HPP


#include "FockSpace/FockSpaceType.hpp"
#include "math/Matrix.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "EvaluationContainer.hpp"

#include <iostream>
#include <memory>



namespace GQCP {


/**
 *  A base class for the representation of a Fock space
 */
class BaseFockSpace {
protected:
    size_t K;  // number of spatial orbitals
    size_t dim;  // dimension of the Fock space


    // PROTECTED CONSTRUCTORS
    BaseFockSpace() = default;
    /**
     *  @param K        the number of orbitals
     *  @param dim      the dimension of the Fock space
     */
    BaseFockSpace(size_t K, size_t dim);

public:
    // NAMED CONSTRUCTORS
    /**
     *  Clones a derived BaseFockSpace instance to the heap memory
     *
     *  @param fock_space     reference to a derived BaseFockSpace instance to be cloned.
     *
     *  @return a shared pointer owning the heap-cloned Fock space
     */
    static std::shared_ptr<BaseFockSpace> CloneToHeap(const BaseFockSpace& fock_space);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // GETTERS
    size_t get_dimension() const { return dim; }
    size_t get_K() const { return K; }
    virtual FockSpaceType get_type() const = 0;


    // PUBLIC METHODS
    /**
     *  @return the coefficient vector for the Hartree-Fock wave function (i.e. the 'first' ONV/Slater determinant)
     */
    VectorX<double> HartreeFockExpansion() const;

    /**
     *  @return a random normalized coefficient vector, with coefficients uniformly distributed in [-1, 1]
     */
    VectorX<double> randomExpansion() const;

    /**
     *  @return a constant normalized coefficients vector (i.e. all the coefficients are equal)
     */
    VectorX<double> constantExpansion() const;


    // Virtual
    /**
    *  Evaluate the operator in a dense matrix
    *
    *  @param one_op               the one-electron operator to be evaluated in the Fock space
    *  @param diagonal_values      bool to indicate if diagonal values will be calculated
    *
    *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
    */
    virtual SquareMatrix<double> evaluateOperatorDense(const OneElectronOperator<double>& one_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const OneElectronOperator<double>& one_op,
                                                               bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const TwoElectronOperator<double>& two_op, bool diagonal_values) const = 0;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const TwoElectronOperator<double>& two_op,
                                                               bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
     */
    virtual SquareMatrix<double> evaluateOperatorDense(const HamiltonianParameters<double>& ham_par,
                                                       bool diagonal_values) const = 0;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    virtual Eigen::SparseMatrix<double> evaluateOperatorSparse(const HamiltonianParameters<double>& ham_par,
                                                               bool diagonal_values) const = 0;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const OneElectronOperator<double>& one_op) const = 0;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const TwoElectronOperator<double>& two_op) const = 0;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
     */
    virtual VectorX<double> evaluateOperatorDiagonal(const HamiltonianParameters<double>& ham_par) const = 0;

};


}  // namespace GQCP


#endif  // GQCP_BASEFOCKSPACE_HPP
