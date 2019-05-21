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
#ifndef GQCP_PRODUCTFOCKSPACE_HPP
#define GQCP_PRODUCTFOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


/**
 *  A class that represents the product of two full Fock spaces (alpha and beta).
 */
class ProductFockSpace: public BaseFockSpace {
private:
    FockSpace fock_space_alpha;
    FockSpace fock_space_beta;

    std::vector<Eigen::SparseMatrix<double>> alpha_couplings;

    // PRIVATE METHODS
    OneElectronOperator<double> oneElectronPartition(size_t p, size_t q, const TwoElectronOperator<double>& two_op) const;

public:
    // CONSTRUCTORS
    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     */
    ProductFockSpace(size_t K, size_t N_alpha, size_t N_beta);


    // DESTRUCTORS
    ~ProductFockSpace() override = default;


    // GETTERS
    size_t get_N_alpha() const { return this->fock_space_alpha.get_N(); }
    size_t get_N_beta() const { return this->fock_space_beta.get_N(); }
    const FockSpace& get_fock_space_alpha() const { return this->fock_space_alpha; }
    const FockSpace& get_fock_space_beta() const { return this->fock_space_beta; }
    FockSpaceType get_type() const override { return FockSpaceType::ProductFockSpace; }


    // STATIC PUBLIC METHODS
    /**
     *  @param K            the number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     *
     *  @return the dimension of the product Fock space
     */
    static size_t calculateDimension(size_t K, size_t N_alpha, size_t N_beta);


    // PUBLIC METHODS
    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> EvaluateOperatorDense(const OneElectronOperator<double>& one_op, bool diagonal_values = true) const;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> EvaluateOperatorSparse(const OneElectronOperator<double>& one_op, bool diagonal_values = true) const;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> EvaluateOperatorDense(const TwoElectronOperator<double>& two_op, bool diagonal_values = true) const;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> EvaluateOperatorSparse(const TwoElectronOperator<double>& two_op, bool diagonal_values = true) const;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> EvaluateOperatorDense(const HamiltonianParameters<double>& ham_par, bool diagonal_values = true) const;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> EvaluateOperatorSparse(const HamiltonianParameters<double>& ham_par, bool diagonal_values = true) const;

};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
