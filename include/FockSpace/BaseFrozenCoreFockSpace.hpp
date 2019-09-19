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


#include "Basis/SingleParticleBasis.hpp"
#include "FockSpace/BaseFockSpace.hpp"


namespace GQCP {

/**
 *  Struct to accommodate the "freezeOperator(ScalarSQTwoElectronOperator<double>)" method, as it returns a one- and two-electron operator
 */
struct FrozenOperators {
    ScalarSQOneElectronOperator<double> one_op;
    ScalarSQTwoElectronOperator<double> two_op;
};


/**
 *  A base class that for a "frozen" Fock space: this is a subspace in which the first X orbitals are always occupied
 */
class BaseFrozenCoreFockSpace: public BaseFockSpace {
protected:
    size_t X;  // number of frozen orbitals/electrons
    std::shared_ptr<BaseFockSpace> active_fock_space;  // active (non-frozen) Fock space containing only the active electrons (N-X) and orbitals (K-X)

public:
    // CONSTRUCTORS
    /**
     *  @param fock_space                    shared pointer to active (non-frozen core) Fock space
     *  @param X                             the number of frozen orbitals
     */
    BaseFrozenCoreFockSpace(std::shared_ptr<BaseFockSpace> fock_space, size_t X);


    // PUBLIC OVERRIDEN METHODS
    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;
    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;
    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;
    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;
    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const override;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param sq_hamiltonian              the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
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
     *  @return the operator diagonal from strictly evaluating the frozen orbitals in the Fock space
     */
    static VectorX<double> frozenCoreDiagonal(const ScalarSQOneElectronOperator<double>& one_op, size_t X, size_t dimension);

    /**
     *  @param two_op       the two-electron operator in an orthonormal orbital basis
     *  @param X            the number of frozen orbitals
     *
     *  @return the operator diagonal from strictly evaluating the frozen orbitals in the Fock space
     */
    static VectorX<double> frozenCoreDiagonal(const ScalarSQTwoElectronOperator<double>& two_op, size_t X, size_t dimension);

    /**
     *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
     *  @param X                    the number of frozen orbitals
     *  @param dimension            the dimension of the diagonal
     *
     *  @return the Hamiltonian diagonal from strictly evaluating the frozen orbitals in a (any) Fock space
     */
    static VectorX<double> frozenCoreDiagonal(const SQHamiltonian<double>& sq_hamiltonian, size_t X, size_t dimension);
};


}  // namespace GQCP
