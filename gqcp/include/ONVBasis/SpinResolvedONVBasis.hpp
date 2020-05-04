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


#include "ONVBasis/BaseONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A full spin-resolved spin-resolved ONV basis.
 */
class SpinResolvedONVBasis: public BaseONVBasis {
private:
    SpinUnresolvedONVBasis fock_space_alpha;
    SpinUnresolvedONVBasis fock_space_beta;

    std::vector<Eigen::SparseMatrix<double>> alpha_couplings;


public:
    // CONSTRUCTORS

    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     */
    SpinResolvedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta);


    // DESTRUCTORS

    ~SpinResolvedONVBasis() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->fock_space_alpha.get_N(); }
    size_t get_N_beta() const { return this->fock_space_beta.get_N(); }
    const SpinUnresolvedONVBasis& get_fock_space_alpha() const { return this->fock_space_alpha; }
    const SpinUnresolvedONVBasis& get_fock_space_beta() const { return this->fock_space_beta; }
    ONVBasisType get_type() const override { return ONVBasisType::SpinResolvedONVBasis; }
    const std::vector<Eigen::SparseMatrix<double>>& get_alpha_couplings() const { return alpha_couplings; }


    // STATIC PUBLIC METHODS

    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     *
     *  @return the dimension of this ONV basis
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);


    // PUBLIC METHODS

    /**
     *  @return the dimension of this ONV basis
     */
    size_t dimension() const;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of this basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of this ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;

    /**
     *  Evaluate the unrestricted Hamiltonian in a dense matrix
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param diagonal_values                bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of this ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const override;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param sq_hamiltonian              the Hamiltonian expressed in an orthonormal basis
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const override;

    /**
     *  Evaluate the diagonal of the unrestricted Hamiltonian
     *
     *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const;

    /**
     *  Evaluate a one electron operator in a matrix vector product
     *
     *  @param one_op                       the one electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in this ONV basis
     *
     *  @return the one electron operator's matrix vector product in a vector with the dimensions of this ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate a two electron operator in a matrix vector product
     *
     *  @param two_op                       the two electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in this ONV basis
     *
     *  @return the two electron operator's matrix vector product in a vector with the dimensions of this ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the Hamiltonian in a matrix vector product
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in this ONV basis
     *
     *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of this ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the unrestricted Hamiltonian in a matrix vector product
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param x                              the vector upon which the evaluation acts 
     *  @param diagonal                       the diagonal evaluated in this ONV basis
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in this ONV basis
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of this ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of this ONV basis
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const override;

    /**
     *  Auxiliary method in order to calculate "theta(pq)",
     *  it returns a partition of a two-electron operator as one-electron operator
     *  where A (i,j) = T (p, q, i, j).
     *
     *  @param p            first fixed index of the two-electron operator
     *  @param q            second fixed index of the two-electron operator
     *  @param two_op       the two-electron operator
     *
     *  @return a one-electron operator containing a partition of the two-electron operator
     */
    ScalarSQOneElectronOperator<double> oneElectronPartition(size_t p, size_t q, const ScalarSQTwoElectronOperator<double>& two_op) const;
};


}  // namespace GQCP
