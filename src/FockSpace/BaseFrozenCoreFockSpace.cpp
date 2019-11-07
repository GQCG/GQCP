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
#include "FockSpace/BaseFrozenCoreFockSpace.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space                    shared pointer to active (non-frozen core) Fock space
 *  @param X                             the number of frozen orbitals
 */
BaseFrozenCoreFockSpace::BaseFrozenCoreFockSpace(std::shared_ptr<GQCP::BaseFockSpace> fock_space, size_t X) :
    BaseFockSpace(fock_space->get_K()+X, fock_space->get_dimension()),
    active_fock_space (std::move(fock_space)),
    X (X)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    // Freeze the Hamiltonian
    ScalarSQOneElectronOperator<double> frozen_one_op = BaseFrozenCoreFockSpace::freezeOperator(one_op, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorDense(frozen_one_op, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(one_op, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    // Freeze the operator
    ScalarSQOneElectronOperator<double> frozen_one_op = BaseFrozenCoreFockSpace::freezeOperator(one_op, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorSparse(frozen_one_op, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(one_op, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    // Freeze the operators
    const auto frozen_ops = BaseFrozenCoreFockSpace::freezeOperator(two_op, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorDense(frozen_ops.one_op, diagonal_values);
    evaluation += this->active_fock_space->evaluateOperatorDense(frozen_ops.two_op, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(two_op, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    // Freeze the operators
    const auto frozen_ops = BaseFrozenCoreFockSpace::freezeOperator(two_op, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorSparse(frozen_ops.one_op, diagonal_values);
    evaluation += this->active_fock_space->evaluateOperatorSparse(frozen_ops.two_op, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(two_op, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {
    // Freeze the operators
    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorDense(frozen_ham_par, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> BaseFrozenCoreFockSpace::evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {

    // Freeze the operators
    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, this->X);

    // Evaluate the frozen operator in the active space
    auto evaluation = this->active_fock_space->evaluateOperatorSparse(frozen_ham_par, diagonal_values);

    if (diagonal_values) {
        // Diagonal correction
        const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian, this->X, this->active_fock_space->get_dimension());
        evaluation += frozen_core_diagonal.asDiagonal();
    }

    return evaluation;
}


/**
 *  Evaluate the diagonal of the operator
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const auto frozen_op = BaseFrozenCoreFockSpace::freezeOperator(one_op, this->X);

    // Calculate diagonal in the active space with the "frozen" Hamiltonian
    const auto diagonal = this->active_fock_space->evaluateOperatorDiagonal(frozen_op);

    // Calculate diagonal for the frozen orbitals
    const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(one_op, this->X, this->active_fock_space->get_dimension());

    return diagonal + frozen_core_diagonal;
};


/**
 *  Evaluate the diagonal of the operator
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto frozen_ops = BaseFrozenCoreFockSpace::freezeOperator(two_op, this->X);

    // Calculate diagonal in the active space with the "frozen" Hamiltonian
    auto diagonal = this->active_fock_space->evaluateOperatorDiagonal(frozen_ops.two_op);
    diagonal += this->active_fock_space->evaluateOperatorDiagonal(frozen_ops.one_op);

    // Calculate diagonal for the frozen orbitals
    const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(two_op, this->X, this->active_fock_space->get_dimension());

    return diagonal + frozen_core_diagonal;
}


/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {

    const auto frozen_ham_par = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian, this->X);

    // Calculate diagonal in the active space with the "frozen" Hamiltonian
    const auto diagonal = this->active_fock_space->evaluateOperatorDiagonal(frozen_ham_par);

    // Calculate diagonal for the frozen orbitals
    const auto frozen_core_diagonal = BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian, this->X, this->active_fock_space->get_dimension());

    return diagonal + frozen_core_diagonal;
}


/**
 *  @param one_op       the one-electron operator in an orthonormal orbital basis
 *  @param X            the number of frozen orbitals
 *
 *  @return 'frozen' one-electron operator which cover evaluations from the active and inactive orbitals
 */
ScalarSQOneElectronOperator<double> BaseFrozenCoreFockSpace::freezeOperator(const ScalarSQOneElectronOperator<double>& one_op, size_t X) {

    size_t K_active = one_op.dimension() - X;
    return ScalarSQOneElectronOperator<double>({one_op.parameters().block(X, X, K_active, K_active)});
}


/**
 *  @param two_op       the two-electron operator in an orthonormal orbital basis
 *  @param X            the number of frozen orbitals
 *
 *  @return 'frozen' two-electron operators as a struct of a one- and two-electron operator which cover evaluations from the active and inactive orbitals
 */
FrozenOperators BaseFrozenCoreFockSpace::freezeOperator(const ScalarSQTwoElectronOperator<double>& two_op, size_t X) {

    size_t K_active = two_op.dimension() - X;
    QCMatrix<double> frozen_one_op_par = QCMatrix<double>::Zero(K_active, K_active);
    const auto& two_op_par = two_op.parameters();
    const auto frozen_two_op_par = QCRankFourTensor<double>::FromBlock(two_op_par, X, X, X, X);

    // Frozen two-electron integrals can be rewritten partially as one electron integrals.
    for (size_t i = 0; i < K_active; i++) {  // iterate over the active orbitals
        size_t q = i + X;  // map active orbitals indexes to total orbital indexes (those including the frozen orbitals)

        for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals
            frozen_one_op_par(i,i) += two_op_par(q,q,l,l) + two_op_par(l,l,q,q) - two_op_par(q,l,l,q)/2 - two_op_par(l,q,q,l)/2;
        }

        for (size_t j = i+1; j < K_active; j++) {  // iterate over the active orbitals
            size_t p = j + X;  // map active orbitals indexes to total orbital indexes (those including the frozen orbitals)

            for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals
                frozen_one_op_par(i,j) += two_op_par(q,p,l,l) + two_op_par(l,l,q,p) - two_op_par(q,l,l,p)/2 - two_op_par(l,p,q,l)/2;
                frozen_one_op_par(j,i) += two_op_par(p,q,l,l) + two_op_par(l,l,p,q) - two_op_par(p,l,l,q)/2 - two_op_par(l,q,p,l)/2;
            }
        }
    }

    return {ScalarSQOneElectronOperator<double>({frozen_one_op_par}), ScalarSQTwoElectronOperator<double>({frozen_two_op_par})};
}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param sp_basis             the single-particle basis
 *  @param X                    the number of frozen orbitals
 *
 *  @return a 'frozen' Hamiltonian which cover two-electron integral evaluations from the active and inactive orbitals
 *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
 */
SQHamiltonian<double> BaseFrozenCoreFockSpace::freezeOperator(const SQHamiltonian<double>& sq_hamiltonian, size_t X) {

    size_t K_active = sq_hamiltonian.dimension() - X;  // number of non-frozen orbitals

    const auto frozen_components_g = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian.twoElectron(), X);

    ScalarSQOneElectronOperator<double> h = BaseFrozenCoreFockSpace::freezeOperator(sq_hamiltonian.core(), X) + frozen_components_g.one_op;  // active
    ScalarSQTwoElectronOperator<double> g = frozen_components_g.two_op;

    return SQHamiltonian<double>(h, g);
}


/**
 *  @param one_op       the one-electron operator in an orthonormal orbital basis
 *  @param X            the number of frozen orbitals
 *
 *  @return the operator diagonal from strictly evaluating the frozen orbitals in the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::frozenCoreDiagonal(const ScalarSQOneElectronOperator<double>& one_op, size_t X, size_t dimension) {

    const auto& one_op_par = one_op.parameters();

    // The diagonal value for the frozen orbitals is the same for each ONV
    double value = 0;
    for (size_t i = 0; i < X; i++) {
        value += 2 * one_op_par(i,i);
    }

    return  VectorX<double>::Constant(dimension, value);
}


/**
 *  @param two_op       the two-electron operator in an orthonormal orbital basis
 *  @param X            the number of frozen orbitals
 *
 *  @return the operator diagonal from strictly evaluating the frozen orbitals in the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::frozenCoreDiagonal(const ScalarSQTwoElectronOperator<double>& two_op, size_t X, size_t dimension) {

    const auto& two_op_par = two_op.parameters();

    // The diagonal value for the frozen orbitals is the same for each ONV
    double value = 0;
    for (size_t i = 0; i < X; i++) {
        value += two_op_par(i,i,i,i);

        for (size_t j = i+1; j < X; j++) {
            value += 2*two_op_par(i,i,j,j) + 2*two_op_par(j,j,i,i) - two_op_par(j,i,i,j) - two_op_par(i,j,j,i);
        }
    }

    return VectorX<double>::Constant(dimension, value);
}


/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param X                    the number of frozen orbitals
 *
 *  @return the Hamiltonian diagonal from strictly evaluating the frozen orbitals in the Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::frozenCoreDiagonal(const SQHamiltonian<double>& sq_hamiltonian, size_t X,  size_t dimension) {
    return BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian.core(), X, dimension) + BaseFrozenCoreFockSpace::frozenCoreDiagonal(sq_hamiltonian.twoElectron(), X, dimension);
}


/*
 *  STATIC UNRESTRICTED
 */ 

/**
 *  @param usq_hamiltonian      the Hamiltonian expressed in an unrestricted orthonormal basis
 *  @param X                    the number of frozen orbitals
 *  @param dimension            the dimension of the diagonal
 *
 *  @return the Hamiltonian diagonal from strictly evaluating the frozen orbitals in a (any) Fock space
 */
VectorX<double> BaseFrozenCoreFockSpace::frozenCoreDiagonal(const USQHamiltonian<double>& usq_hamiltonian, size_t X, size_t dimension) {

    QCMatrix<double> one_op_par_alpha = usq_hamiltonian.alphaHamiltonian().core().parameters();
    QCMatrix<double> one_op_par_beta = usq_hamiltonian.betaHamiltonian().core().parameters();

    const auto& two_op_par_mixed = usq_hamiltonian.twoElectronMixed().parameters();
    const auto& two_op_par_a = usq_hamiltonian.alphaHamiltonian().twoElectron().parameters();
    const auto& two_op_par_b = usq_hamiltonian.betaHamiltonian().twoElectron().parameters();
    
    // The diagonal value for the frozen orbitals is the same for each ONV
    for (size_t i = 0; i < X; i++) {
        value += two_op_par_mixed(i,i,i,i);
        value += one_op_par_alpha(i,i) + one_op_par_beta(i,i);
        for (size_t j = i+1; j < X; j++) {
            value += two_op_par_mixed(i,i,j,j) + two_op_par_mixed(j,j,i,i) + two_op_par_a(i,i,j,j)/2 + two_op_par_a(j,j,i,i)/2 - two_op_par_a(j,i,i,j)/2 - two_op_par_a(i,j,j,i)/2  + two_op_par_b(i,i,j,j)/2 + two_op_par_b(j,j,i,i)/2 - two_op_par_b(j,i,i,j)/2 - two_op_par_b(i,j,j,i)/2;
        }
    }

    return VectorX<double>::Constant(dimension, value);
}


/**
 *  @param usq_hamiltonian      the Hamiltonian expressed in an unrestricted orthonormal basis
 *  @param X                    the number of frozen orbitals
 *
 *  @return a 'frozen' Hamiltonian which cover two-electron integral evaluations from the active and inactive orbitals
 *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
 */
USQHamiltonian<double> BaseFrozenCoreFockSpace::freezeOperator(const USQHamiltonian<double>& usq_hamiltonian, size_t X) {
    
    size_t K_active = usq_hamiltonian.dimension() - X;  // number of non-frozen orbitals

    QCMatrix<double> frozen_one_op_par_alpha = usq_hamiltonian.alphaHamiltonian().core().parameters().block(X, X, K_active, K_active);
    QCMatrix<double> frozen_one_op_par_beta = usq_hamiltonian.betaHamiltonian().core().parameters().block(X, X, K_active, K_active);

    const auto& two_op_par_alpha = usq_hamiltonian.alphaHamiltonian().twoElectron().parameters();
    const auto& two_op_par_beta = usq_hamiltonian.betaHamiltonian().twoElectron().parameters();
    const auto& two_op_par_mixed = usq_hamiltonian.twoElectronMixed().parameters();

    QCRankFourTensor<double> frozen_two_op_par_alpha = QCRankFourTensor<double>::FromBlock(usq_hamiltonian.alphaHamiltonian().twoElectron().parameters(), X, X, X, X);
    const auto frozen_two_op_par_beta = QCRankFourTensor<double>::FromBlock(usq_hamiltonian.betaHamiltonian().twoElectron().parameters(), X, X, X, X);
    const auto frozen_two_op_par_mixed = QCRankFourTensor<double>::FromBlock(usq_hamiltonian.twoElectronMixed().parameters(), X, X, X, X);

    // Frozen two-electron integrals can be rewritten partially as one electron integrals.
    for (size_t i = 0; i < K_active; i++) {  // iterate over the active orbitals
        size_t q = i + X;  // map active orbitals indexes to total orbital indexes (those including the frozen orbitals)

        for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals
            frozen_one_op_par_alpha(i,i) += two_op_par_mixed(q,q,l,l) + two_op_par_alpha(l,l,q,q) - two_op_par_alpha(q,l,l,q)/2 - two_op_par_alpha(l,q,q,l)/2;
            // There's a small catch here, the variable 'two_op_par_mixed' is represented as g_aabb, the formulas dictate we need g_bbaa for the beta component, hence we switched the indices
            frozen_one_op_par_beta(i,i) += two_op_par_mixed(l,l,q,q) + two_op_par_beta(l,l,q,q) - two_op_par_beta(q,l,l,q)/2 - two_op_par_beta(l,q,q,l)/2;
        }

        for (size_t j = i+1; j < K_active; j++) {  // iterate over the active orbitals
            size_t p = j + X;  // map active orbitals indexes to total orbital indexes (those including the frozen orbitals)

            for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals
                frozen_one_op_par_alpha(i,j) += two_op_par_mixed(q,p,l,l) + two_op_par_alpha(l,l,q,p) - two_op_par_alpha(q,l,l,p)/2 - two_op_par_alpha(l,p,q,l)/2;
                frozen_one_op_par_beta(i,j) += two_op_par_mixed(l,l,q,p) + two_op_par_beta(l,l,q,p) - two_op_par_beta(q,l,l,p)/2 - two_op_par_beta(l,p,q,l)/2;

                frozen_one_op_par_alpha(j,i) += two_op_par_mixed(p,q,l,l) + two_op_par_alpha(l,l,p,q) - two_op_par_alpha(p,l,l,q)/2 - two_op_par_alpha(l,q,p,l)/2;
                frozen_one_op_par_beta(j,i) += two_op_par_mixed(l,l,p,q) + two_op_par_beta(l,l,p,q) - two_op_par_beta(p,l,l,q)/2 - two_op_par_beta(l,q,p,l)/2;
            }
        }
    }

    return USQHamiltonian<double>(SQHamiltonian<double>(ScalarSQOneElectronOperator<double>({frozen_one_op_par_alpha}), ScalarSQTwoElectronOperator<double>({frozen_two_op_par_alpha})), SQHamiltonian<double>(ScalarSQOneElectronOperator<double>({frozen_one_op_par_beta}), ScalarSQTwoElectronOperator<double>({frozen_two_op_par_beta})), ScalarSQTwoElectronOperator<double>({frozen_two_op_par_mixed}));
}


}  // namespace GQCP
