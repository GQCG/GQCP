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


#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/Configuration.hpp"
#include "FockSpace/EvaluationIterator.hpp"
#include "FockSpace/FrozenProductFockSpace.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A class that represents a Fock space that is flexible in the number of states that span it
 *
 *  Configurations are represented as a Configuration: a combination of an alpha and a beta ONV
 */
class SelectedFockSpace : public BaseFockSpace {
private:
    size_t N_alpha;  // number of alpha electrons
    size_t N_beta;  // number of beta electrons

    std::vector<Configuration> configurations;

    /**
     *  @param onv1     the alpha ONV as a string representation read from right to left
     *  @param onv2     the beta ONV as a string representation read from right to left
     *
     *  @return the configuration that holds both ONVs
     *
     *  IMPORTANT: only works for up to 64 bits!
     */
    Configuration makeConfiguration(const std::string& onv1, const std::string& onv2) const;

public:
    // CONSTRUCTORS
    SelectedFockSpace() = default;  // need a default constructor

    /**
     *  A constructor with initial Fock space dimension of 0
     *
     *  @param K            the number of orbitals
     *  @param N_alpha      the number of alpha electrons
     *  @param N_beta       the number of beta electrons
     */
    SelectedFockSpace(size_t K, size_t N_alpha, size_t N_beta);

    /**
     *  A constructor that generates the configurations based on the given ProductFockSpace.
     *
     *  @param fock_space       the product Fock space from which the configurations should be generated
     */
    explicit SelectedFockSpace(const ProductFockSpace& fock_space);

    /**
     *  A constructor that generates the configurations based on the given FockSpace.
     *
     *  @param fock_space       the Fock space from which the configurations should be generated
     */
    explicit SelectedFockSpace(const FockSpace& fock_space);

    /**
     *  A constructor that generates the configurations based on the given frozen product Fock space.
     *
     *  @param fock_space       the frozen product Fock space from which the configurations should be generated
     */
    explicit SelectedFockSpace(const FrozenProductFockSpace& fock_space);

    /**
     *  A constructor that generates the configurations based on the given frozen Fock space.
     *
     *  @param fock_space       the frozen Fock space from which the configurations should be generated
     */
    explicit SelectedFockSpace(const FrozenFockSpace& fock_space);


    // GETTERS
    size_t get_N_alpha() const { return this->N_alpha; }
    size_t get_N_beta() const { return this->N_beta; }
    const Configuration& get_configuration(size_t index) const { return this->configurations[index]; }
    FockSpaceType get_type() const override { return FockSpaceType::SelectedFockSpace; }


    // PUBLIC METHODS
    /**
     *  Make a configuration (see makeConfiguration()) and add it to this Fock space
     *
     *  @param onv1     the alpha ONV as a string representation read from right to left
     *  @param onv2     the beta ONV as a string representation read from right to left
     */
    void addConfiguration(const std::string& onv1, const std::string& onv2);

    /**
     *  Make configurations (see makeConfiguration()) and add them to the Fock space
     *
     *  @param onv1s     the alpha ONVs as string representations read from right to left
     *  @param onv2s     the beta ONVs as string representations read from right to left
     */
    void addConfiguration(const std::vector<std::string>& onv1s, const std::vector<std::string>& onv2s);

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

    /**
     *  Evaluate a one electron operator in a matrix vector product
     *
     *  @param one_op                       the one electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the Fock space
     *
     *  @return the one electron operator's matrix vector product in a vector with the dimensions of the Fock space
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate a two electron operator in a matrix vector product
     *
     *  @param two_op                       the two electron operator expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the Fock space
     *
     *  @return the two electron operator's matrix vector product in a vector with the dimensions of the Fock space
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the Hamiltonian in a matrix vector product
     *
     *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the Fock space
     *
     *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of the Fock space
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;


    // UNRESTRICTED METHODS
    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const;

    /**
     *  Evaluate the Hamiltonian in a matrix vector product
     *
     *  @param usq_hamiltonian              the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param x                            the vector upon which the evaluation acts 
     *  @param diagonal                     the diagonal evaluated in the Fock space
     *
     *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of the Fock space
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  Evaluate the Hamiltonian in a sparse matrix
     *
     *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param diagonal_values          bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
     */
    Eigen::SparseMatrix<double> evaluateOperatorSparse(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const;


    
    /**
     *  PUBLIC TEMPLATED METHODS
     */

    /**
     *  Evaluate the operator in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void EvaluateOperator(const ScalarSQOneElectronOperator<double>& one_op, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {
        // Calling the unrestricted universal method, with identical alpha and beta components does not affect the performance, hence we avoid duplicated code for the restricted part.
        EvaluateOperator(one_op, one_op, evaluation_iterator, diagonal_values);
    }


    /**
     *  Evaluate the operator in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param usq_hamiltonian               the Hamiltonian expressed in an unrestricted orthonormal basis
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void EvaluateOperator(const USQHamiltonian<double>& usq_hamiltonian, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {

        if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
            throw std::invalid_argument("SelectedFockSpace::EvaluateOperator(USQHamiltonian<double>, EvaluationIterator&, bool): Different spinor dimensions of spin components are currently not supported.");
        }
    
        EvaluateOperator(usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA).core(), usq_hamiltonian.spinHamiltonian(SpinComponent::BETA).core(), usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA).twoElectron(), usq_hamiltonian.spinHamiltonian(SpinComponent::BETA).twoElectron(), usq_hamiltonian.twoElectronMixed(), evaluation_iterator, diagonal_values);
    }


    /**
     *  Evaluate the operator in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op_alpha                  the alpha component of a one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param one_op_beta                   the beta component of a one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void EvaluateOperator(const ScalarSQOneElectronOperator<double>& one_op_alpha, const ScalarSQOneElectronOperator<double>& one_op_beta, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {

        const size_t dim = this->get_dimension();
        const auto& h_a = one_op_alpha.parameters();
        const auto& h_b = one_op_beta.parameters();

        for (;!evaluation_iterator.is_finished(); evaluation_iterator.increment()) {  // loop over all addresses (1)
            Configuration configuration_I = this->get_configuration(evaluation_iterator.index);
            ONV alpha_I = configuration_I.onv_alpha;
            ONV beta_I = configuration_I.onv_beta;

            if (diagonal_values) {
                for (size_t p = 0; p < K; p++) {
                    if (alpha_I.isOccupied(p)) {
                        evaluation_iterator.addRowwise(evaluation_iterator.index, h_a(p, p));
                    }

                    if (beta_I.isOccupied(p)) {
                        evaluation_iterator.addRowwise(evaluation_iterator.index, h_b(p,p));
                    }
                }  // loop over q
            }

            // Calculate the off-diagonal elements, by going over all other ONVs
            for (size_t J = evaluation_iterator.index+1; J < dim; J++) {

                Configuration configuration_J = this->get_configuration(J);
                ONV alpha_J = configuration_J.onv_alpha;
                ONV beta_J = configuration_J.onv_beta;

                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(
                            alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(
                            alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);

                    double value = h_a(p, q);

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);
                }

                // 0 electron excitations in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {


                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);

                    double value = h_b(p,q);

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);
                }
            }  // loop over addresses J > I
        }  // loop over addresses I
    }


    /**
     *  Evaluate the operator in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
   template <typename _Matrix>
    void EvaluateOperator(const ScalarSQTwoElectronOperator<double>& two_op, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {
        // Calling this combined method for both the one- and two-electron operator does not affect the performance, hence we avoid writting more code by plugging a zero one-electron operator in the combined method.
        EvaluateOperator(ScalarSQOneElectronOperator<double>(this->K), ScalarSQOneElectronOperator<double>(this->K), two_op, two_op, two_op, evaluation_iterator, diagonal_values);
    }


    /**
     *  Evaluate the operator in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op                        the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param two_op                        the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param evaluation_iterator           evaluation iterator to which the evaluations are added
     *  @param diagonal_values               bool to indicate if diagonal values will be calculated
     */
   template <typename _Matrix>
    void EvaluateOperator(const ScalarSQOneElectronOperator<double>& one_op, const ScalarSQTwoElectronOperator<double>& two_op, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {
        // Calling the unrestricted universal method, with identical alpha, beta and mixed components does not affect the performance, hence we avoid duplicated code for the restricted part.
        EvaluateOperator(one_op, one_op, two_op, two_op, two_op, evaluation_iterator, diagonal_values);
    }


    /**
     *  Evaluate the operators in a given evaluation iterator in the Fock space
     *
     *  @tparam Matrix                       the type of matrix used to store the evaluations
     *
     *  @param one_op_alpha                     the alpha component of a one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param one_op_beta                      the beta component of a one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param two_op_alpha                     the alpha component of a two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param two_op_beta                      the beta component of a two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param two_op_mixed                     the alpha-beta component of a two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
     *  @param evaluation_iterator              evaluation iterator to which the evaluations are added
     *  @param diagonal_values                  bool to indicate if diagonal values will be calculated
     */
    template <typename _Matrix>
    void EvaluateOperator(const ScalarSQOneElectronOperator<double>& one_op_alpha, const ScalarSQOneElectronOperator<double>& one_op_beta, const ScalarSQTwoElectronOperator<double>& two_op_alpha, const ScalarSQTwoElectronOperator<double>& two_op_beta, const ScalarSQTwoElectronOperator<double>& two_op_mixed, EvaluationIterator<_Matrix>& evaluation_iterator, bool diagonal_values) const {

        const size_t dim = this->get_dimension();
        const size_t K = this->get_K();

        const auto& h_a = one_op_alpha.parameters();
        const auto& g_a = two_op_alpha.parameters();
        const auto& h_b = one_op_beta.parameters();
        const auto& g_b = two_op_beta.parameters();

        // Only g_ab is stored, for integrals derived from g_ba we reverse the indices as follows : g_ab(pqrs) = g_ba(rspq)
        const auto& g_ab = two_op_mixed.parameters();

        for ( ;!evaluation_iterator.is_finished(); evaluation_iterator.increment()) {  // loop over all addresses (1)
            Configuration configuration_I = this->get_configuration(evaluation_iterator.index);
            ONV alpha_I = configuration_I.onv_alpha;
            ONV beta_I = configuration_I.onv_beta;

            if (diagonal_values) {
                for (size_t p = 0; p < K; p++) {
                    if (alpha_I.isOccupied(p)) {
                        evaluation_iterator.addRowwise(evaluation_iterator.index, h_a(p,p));
                        for (size_t q = 0; q < K; q++) {

                            if (p != q) {  // can't create/annihilate the same orbital twice
                                if (alpha_I.isOccupied(q)) {
                                    evaluation_iterator.addRowwise(evaluation_iterator.index,  0.5 * g_a(p,p,q,q));
                                    evaluation_iterator.addRowwise(evaluation_iterator.index, -0.5 * g_a(p,q,q,p));
                                }
                            }

                            if (beta_I.isOccupied(q)) {
                                evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * g_ab(p,p,q,q));
                            }
                        }  // loop over q
                    }

                    if (beta_I.isOccupied(p)) {
                        evaluation_iterator.addRowwise(evaluation_iterator.index, h_b(p,p));
                        for (size_t q = 0; q < K; q++) {

                            if (p != q) {  // can't create/annihilate the same orbital twice
                                if (beta_I.isOccupied(q)) {
                                    evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * g_b(p,p,q,q));
                                    evaluation_iterator.addRowwise(evaluation_iterator.index, -0.5 * g_b(p,q,q,p));
                                }
                            }

                            if (alpha_I.isOccupied(q)) {
                                evaluation_iterator.addRowwise(evaluation_iterator.index, 0.5 * g_ab(q,q,p,p));  // g_ab(pqrs) = g_ba(rspq)
                            }
                        }  // loop over q
                    }
                }  // loop over q
            }

            // Calculate the off-diagonal elements, by going over all other ONVs
            for (size_t J = evaluation_iterator.index+1; J < dim; J++) {

                Configuration configuration_J = this->get_configuration(J);
                ONV alpha_J = configuration_J.onv_alpha;
                ONV beta_J = configuration_J.onv_beta;

                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q);

                    double value = h_a(p,q);

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);

                    for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                        if (alpha_I.isOccupied(r) && alpha_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital

                                double value = 0.5 * (g_a(p,q,r,r)
                                                      - g_a(r,q,p,r)
                                                      - g_a(p,r,r,q)
                                                      + g_a(r,r,p,q));

                                evaluation_iterator.addColumnwise(J, sign * value);
                                evaluation_iterator.addRowwise(J, sign * value);
                            }
                        }

                        if (beta_I.isOccupied(r)) {  // beta_I == beta_J from the previous if-branch

                            double value = 0.5 * 2 * g_ab(p,q,r,r);  // g_ab(pqrs) = g_ba(rspq)

                            evaluation_iterator.addColumnwise(J, sign * value);
                            evaluation_iterator.addRowwise(J, sign * value);
                        }
                    }
                }

                // 0 electron excitations in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {


                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    // Calculate the total sign
                    int sign = beta_I.operatorPhaseFactor(p) * beta_J.operatorPhaseFactor(q);

                    double value = h_b(p,q);

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);

                    for (size_t r = 0; r < K; r++) {  // r loops over spatial orbitals

                        if (beta_I.isOccupied(r) && beta_J.isOccupied(r)) {  // r must be occupied on the left and on the right
                            if ((p != r) && (q != r)) {  // can't create or annihilate the same orbital
                                double value = 0.5 * (g_b(p,q,r,r)
                                                      - g_b(r,q,p,r)
                                                      - g_b(p,r,r,q)
                                                      + g_b(r,r,p,q));

                                evaluation_iterator.addColumnwise(J, sign * value);
                                evaluation_iterator.addRowwise(J, sign * value);
                            }
                        }

                        if (alpha_I.isOccupied(r)) {  // alpha_I == alpha_J from the previous if-branch

                            double value = 0.5 * 2 * g_ab(r,r,p,q);  // g_ab(pqrs) = g_ba(rspq)  

                            evaluation_iterator.addColumnwise(J, sign * value);
                            evaluation_iterator.addRowwise(J, sign * value);
                        }
                    }
                }

                // 1 electron excitation in alpha, 1 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 2) && (beta_I.countNumberOfDifferences(beta_J) == 2)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    size_t p = alpha_I.findDifferentOccupations(alpha_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t q = alpha_J.findDifferentOccupations(alpha_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    size_t r = beta_I.findDifferentOccupations(beta_J)[0];  // we're sure that there is only 1 element in the std::vector<size_t>
                    size_t s = beta_J.findDifferentOccupations(beta_I)[0];  // we're sure that there is only 1 element in the std::vector<size_t>

                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_J.operatorPhaseFactor(q) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(s);
                    double value = 0.5 * 2 * g_ab(p,q,r,s);  // g_ab(pqrs) = g_ba(rspq)

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);
                }

                // 2 electron excitations in alpha, 0 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 4) && (beta_I.countNumberOfDifferences(beta_J) == 0)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = alpha_I.findDifferentOccupations(alpha_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = alpha_J.findDifferentOccupations(alpha_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];

                    int sign = alpha_I.operatorPhaseFactor(p) * alpha_I.operatorPhaseFactor(r) * alpha_J.operatorPhaseFactor(q) * alpha_J.operatorPhaseFactor(s);

                    double value = 0.5 * (g_a(p,q,r,s)
                                          - g_a(p,s,r,q)
                                          - g_a(r,q,p,s)
                                          + g_a(r,s,p,q));

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);
                }

                // 0 electron excitations in alpha, 2 in beta
                if ((alpha_I.countNumberOfDifferences(alpha_J) == 0) && (beta_I.countNumberOfDifferences(beta_J) == 4)) {

                    // Find the orbitals that are occupied in one string, and aren't in the other
                    std::vector<size_t> occupied_indices_I = beta_I.findDifferentOccupations(beta_J);  // we're sure this has two elements
                    size_t p = occupied_indices_I[0];
                    size_t r = occupied_indices_I[1];

                    std::vector<size_t> occupied_indices_J = beta_J.findDifferentOccupations(beta_I);  // we're sure this has two elements
                    size_t q = occupied_indices_J[0];
                    size_t s = occupied_indices_J[1];

                    int sign = beta_I.operatorPhaseFactor(p) * beta_I.operatorPhaseFactor(r) * beta_J.operatorPhaseFactor(q) * beta_J.operatorPhaseFactor(s);

                    double value = 0.5 * (g_b(p,q,r,s)
                                          - g_b(p,s,r,q)
                                          - g_b(r,q,p,s)
                                          + g_b(r,s,p,q));

                    evaluation_iterator.addColumnwise(J, sign * value);
                    evaluation_iterator.addRowwise(J, sign * value);
                }
            }  // loop over addresses J > I
        }  // loop over addresses I
    }
};


}  // namespace GQCP
