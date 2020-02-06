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


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "ONVBasis/BaseONVBasis.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A small struct to accommodate the return type of "freezeOperator(const ScalarSQTwoElectronOperator<double>&)", as it should both return a one- and two-electron operator
 */
struct FrozenOperators {
    ScalarSQOneElectronOperator<double> one_op;
    ScalarSQTwoElectronOperator<double> two_op;
};


/**
 *  A base class for a "frozen" ONV basis: this represents an ONV basis in which the first X occupation numbers are always 1
 */
class BaseFrozenCoreONVBasis: public BaseONVBasis {
protected:
    size_t X;  // the number of frozen orbitals/electrons
    std::shared_ptr<BaseONVBasis> active_fock_space;  // active (non-frozen) ONV basis containing only the active electrons (N-X) and orbitals (K-X)

public:

    // CONSTRUCTORS

    /**
     *  @param fock_space                    shared pointer to active (non-frozen core) ONV basis
     *  @param X                             the number of frozen orbitals
     */
    BaseFrozenCoreONVBasis(std::shared_ptr<BaseONVBasis> fock_space, size_t X);


    // PUBLIC OVERRIDEN METHODS

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const override;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param sq_hamiltonian              the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;


    // STATIC PUBLIC METHODS

    /**
     *  @param one_op       the one-electron operator in an orthonormal orbital basis
     *  @param X            the number of frozen orbitals
     *
     *  @return 'frozen' one-electron operator which cover evaluations from the active and inactive orbitals
     */
    static ScalarSQOneElectronOperator<double> freezeOperator(const ScalarSQOneElectronOperator<double>& one_op, size_t X);

    /**
     *  @param two_op       the two-electron operator in an orthonormal orbital basis
     *  @param X            the number of frozen orbitals
     *
     *  @return 'frozen' two-electron operators as a struct of a one- and two-electron operator which cover evaluations from the active and inactive orbitals
     */
    static FrozenOperators freezeOperator(const ScalarSQTwoElectronOperator<double>& two_op, size_t X);


    // STATIC PUBLIC METHODS
    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param X                    the number of frozen orbitals
     *
     *  @return a 'frozen' Hamiltonian which cover two-electron integral evaluations from the active and inactive orbitals
     *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
     */
    static SQHamiltonian<double> freezeOperator(const SQHamiltonian<double>& sq_hamiltonian, size_t X);

    /**
     *  @param one_op       the one-electron operator in an orthonormal orbital basis
     *  @param X            the number of frozen orbitals
     *
     *  @return the operator diagonal from strictly evaluating the frozen orbitals in the ONV basis
     */
    static VectorX<double> frozenCoreDiagonal(const ScalarSQOneElectronOperator<double>& one_op, size_t X, size_t dimension);

    /**
     *  @param two_op       the two-electron operator in an orthonormal orbital basis
     *  @param X            the number of frozen orbitals
     *
     *  @return the operator diagonal from strictly evaluating the frozen orbitals in the ONV basis
     */
    static VectorX<double> frozenCoreDiagonal(const ScalarSQTwoElectronOperator<double>& two_op, size_t X, size_t dimension);

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param X                    the number of frozen orbitals
     *  @param dimension            the dimension of the diagonal
     *
     *  @return the Hamiltonian diagonal from strictly evaluating the frozen orbitals in a (any) ONV basis
     */
    static VectorX<double> frozenCoreDiagonal(const SQHamiltonian<double>& sq_hamiltonian, size_t X, size_t dimension);


    // STATIC UNRESTRICTED
    /**
     *  @param usq_hamiltonian      the Hamiltonian expressed in an unrestricted orthonormal basis
     *  @param X                    the number of frozen orbitals
     *  @param dimension            the dimension of the diagonal
     *
     *  @return the Hamiltonian diagonal from strictly evaluating the frozen orbitals in a (any) ONV basis
     */
    static VectorX<double> frozenCoreDiagonal(const USQHamiltonian<double>& usq_hamiltonian, size_t X, size_t dimension);

    /**
     *  @param usq_hamiltonian      the Hamiltonian expressed in an unrestricted orthonormal basis
     *  @param X                    the number of frozen orbitals
     *
     *  @return a 'frozen' Hamiltonian which cover two-electron integral evaluations from the active and inactive orbitals
     *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
     */
    static USQHamiltonian<double> freezeOperator(const USQHamiltonian<double>& usq_hamiltonian, size_t X);
};


}  // namespace GQCP
