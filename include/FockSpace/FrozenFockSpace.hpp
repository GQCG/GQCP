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
#ifndef GQCP_FROZENFOCKSPACE_HPP
#define GQCP_FROZENFOCKSPACE_HPP

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/FockPermutator.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FrozenCoreFockSpace.hpp"


namespace GQCP {


/**
 *  A class that represents a frozen Fock space: this is a subspace of the N-electron Fock space in which the first X orbitals are always occupied
 */
class FrozenFockSpace: public FrozenCoreFockSpace, public FockPermutator<FrozenFockSpace> {
protected:
    size_t X;  // number of frozen orbitals/electrons
    FockSpace active_fock_space;  // active (non-frozen) Fock space containing only the active electrons (N-X) and orbitals (K-X)

public:
    // CONSTRUCTORS
    /**
     *  @param K        the total number of orbitals
     *  @param N        the total number of electrons
     *  @param X        the number of frozen orbitals and electrons
     */
    FrozenFockSpace(size_t K, size_t N, size_t X);

    /**
     *  @param fock_space       (to be frozen) full Fock space
     *  @param X                the number of frozen orbitals and electrons
     */
    FrozenFockSpace(const FockSpace& fock_space, size_t X);


    // DESTRUCTOR
    ~FrozenFockSpace() override = default;


    // GETTERS
    size_t get_number_of_frozen_orbitals() const { return this->X; }
    const FockSpace& get_active_fock_space() const { return this->active_fock_space; }
    FockSpaceType get_type() const override { return FockSpaceType::FrozenFockSpace; }


    // OVERRIDEN PUBLIC METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation in the frozen Fock space
     *
     *      Examples (first orbital is frozen):
     *          0101 -> 1001
     *         01101 -> 10011
     */
    size_t ulongNextPermutation(size_t representation) const override;

    /**
     *  @param representation      a representation of an ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t calculateRepresentation(size_t address) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const ONV& onv) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const ONV& onv) const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the Fock space
     */
    size_t countTotalOneElectronCouplings() const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the Fock space
     */
    size_t countTotalTwoElectronCouplings() const override;



    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const OneElectronOperator<double>& one_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const OneElectronOperator<double>& one_op,
                                                       bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a dense matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const TwoElectronOperator<double>& two_op, bool diagonal_values) const override;

    /**
     *  Evaluate the operator in a sparse matrix
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const TwoElectronOperator<double>& two_op,
                                                       bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const HamiltonianParameters<double>& ham_par,
                                               bool diagonal_values) const override;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *  @param diagonal_values      bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const HamiltonianParameters<double>& ham_par,
                                                       bool diagonal_values) const override;


    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param one_op               the one-electron operator to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const OneElectronOperator<double>& one_op) const override;

    /**
     *  Evaluate the diagonal of the operator
     *
     *  @param two_op               the two-electron operator to be evaluated in the Fock space
     *
     *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const TwoElectronOperator<double>& two_op) const override;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param ham_par              HamiltonianParameters to be evaluated in the Fock space
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const HamiltonianParameters<double>& ham_par) const override;




    // PUBLIC METHODS
    /**
     *  If we have
     *      FrozenFockSpace fock_space;
     *
     *  This makes sure that we can call
     *      fock_space.getAddress(onv);
     *  instead of the syntax
     *      fock_space.FockPermutator<FrozenFockSpace>::getAddress(onv);
     */
    using FockPermutator<FrozenFockSpace>::getAddress;
};


}  // namespace GQCP


#endif  // GQCP_FROZENFOCKSPACE_HPP
