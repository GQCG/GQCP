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


#include "Basis/BiorthogonalBasis/ULowdinPairingBasis.hpp"
#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "QuantumChemical/Spin.hpp"
#include "Utilities/CRTP.hpp"


namespace GQCP {


/*
 *  MARK: UNonOrthogonalStateBasis
 */

/**
 *  A basis formed by any number of non-orthogonal states, in the form of `UTransformation`s.
 *
 *  @tparam _ExpansionScalar        The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
 */
template <typename _Scalar>
class UNonOrthogonalStateBasis:
    public CRTP<UNonOrthogonalStateBasis<_Scalar>>,
    public BasisTransformable<UNonOrthogonalStateBasis<_Scalar>>,
    public JacobiRotatable<UNonOrthogonalStateBasis<_Scalar>> {

public:
    // The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
    using Scalar = _Scalar;

    // The type of matrix associated with this kind of NonOrthogonalStateBasis.
    using Matrix = SquareMatrix<Scalar>;

    // The type of transformation that is naturally related to a `NonOrthogonalStateBasis`.
    using Transformation = UTransformation<Scalar>;

    // The type of non-orthogonal state basis this is.
    using Self = UNonOrthogonalStateBasis<Scalar>;

    // The biorthogonal basis related to the basis states in this type of non-orthogonal basis.
    using BiorthogonalBasis = ULowdinPairingBasis<Scalar>;

    // The second-quantized representation of any one-electron operator operator related to the `UNonOrthogonalStateBasis`.
    using SQOverlapOperator = ScalarUSQOneElectronOperator<Scalar>;

    // The second-quantized representation of any one-electron operator operator related to the `UNonOrthogonalStateBasis`.
    using OneElectronOperator = ScalarUSQOneElectronOperator<Scalar>;

    // The second-quantized representation of any two-electron operator operator related to the `UNonOrthogonalStateBasis`.
    using TwoElectronOperator = ScalarUSQTwoElectronOperator<Scalar>;

    // The second-quantized representation of the Hamiltonian that can be evaluated in this basis.
    using Hamiltonian = USQHamiltonian<Scalar>;

    // The type of Jacobi rotation that is naturally related to the derived non orthogonal-basis.
    using JacobiRotationType = typename JacobiRotatableTraits<Self>::JacobiRotationType;

    // The vector containing the basis state of the associated type of transformations.
    using States = std::vector<Transformation>;


protected:
    // The vector containing the non-orthogonal basis states.
    States basis_states;

    // The overlap operator in AO basis, constructed from the spinor/spin-orbital basis.
    SQOverlapOperator overlap_operator_AO;

    // The total number of occupied alpha orbitals.
    size_t N_a;

    // The total number of occupied beta orbitals.
    size_t N_b;

    // The threshold used to determine zero values.
    double zero_threshold;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `UNonOrthogonalStateBasis` from any number of non orthogonal states.
     *
     *  @param basis_state_vector                        The vector containing the non-orthogonal basis states.
     *  @param S_AO                                      The overlap operator in AO basis, constructed from the spinor/spin-orbital used to calculate the non-orthogonal states.
     *  @param number_of_occupied_alpha_orbitals         The total number of occupied orbitals in the system.
     *  @param number_of_occupied_beta_orbitals          The total number of occupied orbitals in the system.
     *  @param threshold                                 The threshold at which a value is verified to be zero or not. The default is 1e-8.
     */
    UNonOrthogonalStateBasis<Scalar>(const States& basis_state_vector, const SQOverlapOperator& S_AO, const size_t number_of_occupied_alpha_orbitals, const size_t number_of_occupied_beta_orbitals, const double threshold = 1e-8) :
        basis_states {basis_state_vector},
        overlap_operator_AO {S_AO},
        N_a {number_of_occupied_alpha_orbitals},
        N_b {number_of_occupied_beta_orbitals},
        zero_threshold {threshold} {

        // The basis states must have the same dimensions.
        for (size_t i = 0; i < basis_state_vector.size(); i++) {
            if (basis_state_vector[0].alpha().dimension() != basis_state_vector[i].alpha().dimension()) {
                throw std::invalid_argument("NonOrthogonalStateBasis<Scalar>(const States& basis_state_vector, const SQOverlapOperator& S_AO, const size_t number_of_occupied_orbitals, const double threshold = 1e-8): The given basis states do not have the same dimensions.");
            }
        }
    }


    /*
     *  MARK: Properties
     */

    /**
     * Return the i'th basis states in the formed non-orthogonal state basis.
     *
     * @param i     The index of the basis state requested.
     *
     * @return The i'th basis state..
     */
    const Transformation& basisState(size_t i) const { return this->basis_states[i]; }


    /**
     * Return the dimension of the basis states in the formed non-orthogonal state basis.
     *
     * @return The dimension of the basis states.
     *
     * @note We return the dimension of the first state, as the constructor checks that all states have the same dimension. We also assume the alpha and beta part have the same dimension.
     */
    const size_t basisStateDimension() const { return this->basis_states[0].alpha().dimension(); }

    /**
     * Return the basis states in the formed non-orthogonal state basis.
     *
     * @return The basis states.
     */
    const States& basisStates() const { return this->basis_states; }

    /**
     * Return the number of basis states in the formed non-orthogonal state basis.
     *
     * @return The number of basis states.
     */
    const size_t numberOfBasisStates() const { return this->basis_states.size(); }

    /**
     * Return the threshold used to compare values to zero associated with this non-orthogonal state basis.
     *
     * @return The threshold at which to evaluate zero values..
     */
    const double& threshold() const { return this->zero_threshold; }


    /**
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     *
     *  @param T            The basis transformation.
     *
     *  @return The basis-transformed object.
     */
    Self transformed(const Transformation& T) const override {

        auto result = this->derived();

        for (size_t i = 0; i < result.numberOfBasisStates(); i++) {
            result.basis_states[i].transform(T);
        }

        return result;
    }

    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<Self>::rotated;


    /**
     *  MARK: Conforming to `JacobiRotatable`.
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     *
     *  @param jacobi_rotation          The Jacobi rotation.
     *
     *  @return The Jacobi-rotated object.
     */
    Self rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto Ja = UTransformationComponent<Scalar>::FromJacobi(jacobi_rotation, this->basisStateDimension());
        const auto Jb = UTransformationComponent<Scalar>::FromJacobi(jacobi_rotation, this->basisStateDimension());
        Transformation J {Ja, Jb};

        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<Self>::rotate;


    /**
     * MARK: Operator Evaluations
     */

    /**
     * Evaluate any Hamiltonian operator in this non-orthogonal state basis.
     *
     * @param hamiltonian                       The Hamiltonian operator to be evaluated.
     * @param nuclear_repulsion_operator        The nuclear repulsion operator associated with the molecule used for these calculations.
     *
     * @return The matrix representation of this Hamiltonian expressed in the non-orthogonal state basis.
     *
     * @note The given Hamiltonian must be expressed in the non-orthogonal AO basis.
     * @note In order to express the Hamiltonian in this non-orthogonal state basis, we must take the nuclear repulsion into account as well (see equation (16) in the 2018 paper (https://doi.org/10.1063/1.4999218) by Olsen et. al. for example). This is why we also give the `nuclearRepulsionOperator` as a parameter.
     */
    Matrix evaluateHamiltonianOperator(const Hamiltonian& hamiltonian, const NuclearRepulsionOperator& nuclear_repulsion_operator) const {

        // We first require the separate one and two electron operators from the Hamiltonian we want to evaluate.
        const auto& h = hamiltonian.core();
        const auto& g = hamiltonian.twoElectron();

        // From the nuclear repulsion operator we need the value.
        const auto nuc_rep = nuclear_repulsion_operator.value();

        // Evaluate the operators in this basis.
        // We also need the evaluated overlap in this basis.
        const auto h_evaluated = this->evaluateOneElectronOperator(h);
        const auto g_evaluated = this->evaluateTwoElectronOperator(g);
        const auto s_evaluated = this->evaluateOverlapOperator();

        // Return the Hamiltonian matrix, represented in this non-orthogonal state basis.
        return h_evaluated + g_evaluated + (nuc_rep * s_evaluated);
    }


    /**
     * Evaluate any scalar one-electron operator in this non-orthogonal state basis.
     *
     * @param f_op      The scalar one-electron operator to be evaluated.
     *
     * @return The matrix representation of this operator expressed in the non-orthogonal state basis.
     *
     * @note The given scalar one-electron operator must be expressed in the non-orthogonal AO basis.
     * @note This implementation is based on equation 59 from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     */
    Matrix evaluateOneElectronOperator(const OneElectronOperator& f_op) const {

        // We start by initializing a zero matrix of the correct dimension.
        // Since the matrix is always squared, we only need one parameter to define its dimensions.
        Matrix evaluated_operator = Matrix::Zero(this->numberOfBasisStates());

        // Secondly, we map the operator parameters to a tensor representation.
        // This will be necessary to perform the correct contractions.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> operator_tensor_map_alpha {f_op.alpha().parameters().data(), f_op.alpha().parameters().rows(), f_op.alpha().parameters().cols()};
        Tensor<Scalar, 2> operator_parameters_tensor_alpha = Tensor<Scalar, 2>(operator_tensor_map_alpha);

        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> operator_tensor_map_beta {f_op.beta().parameters().data(), f_op.beta().parameters().rows(), f_op.beta().parameters().cols()};
        Tensor<Scalar, 2> operator_parameters_tensor_beta = Tensor<Scalar, 2>(operator_tensor_map_beta);

        // Loop over all basis states and calculate each matrix element using the generalized Slater-Condon rules.
        for (size_t i = 0; i < this->numberOfBasisStates(); i++) {
            for (size_t j = 0; j < this->numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N_a, this->N_b};

                // Check the number of zeros in the biorthogonal basis, as this influences how the matrix elements are calculated.
                const auto number_of_zeros = lowdin_pairing_basis.numberOfZeroOverlaps();

                // Use the correct formula, depending on how many zero overlaps there are.
                // If there are zero overlap values, we perform the following calculation.
                if (number_of_zeros == 0) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap and the weighted co-density matrix from the biorthogonal basis.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();
                    const auto weighted_co_density_alpha = lowdin_pairing_basis.weightedCoDensity().alpha();
                    const auto weighted_co_density_beta = lowdin_pairing_basis.weightedCoDensity().beta();

                    // Perform the contraction.
                    Tensor<Scalar, 0> matrix_element_alpha = operator_parameters_tensor_alpha.template einsum<2>("uv, vu ->", weighted_co_density_alpha.matrix());
                    Tensor<Scalar, 0> matrix_element_beta = operator_parameters_tensor_beta.template einsum<2>("uv, vu ->", weighted_co_density_beta.matrix());

                    evaluated_operator(i, j) += reduced_overlap * (matrix_element_alpha(0) + matrix_element_beta(0));
                }
                // If there is one zero overlap value, we perform the following calculation.
                else if (number_of_zeros == 1) {

                    // Determine whether the zero stems from the alpha or the beta component.
                    auto zero_spin = Spin::alpha;
                    if (lowdin_pairing_basis.numberOfZeroOverlaps(Spin::beta) > lowdin_pairing_basis.numberOfZeroOverlaps(Spin::alpha)) {
                        zero_spin = Spin::beta;
                    }

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap(zero_spin);

                    // Next, calculate the co-density matrix of the orbital corresponding to the zero overlap value.
                    // We know there is only one zero overlap orbital, so we acces the first index of the vector.
                    const auto zero_overlap_index = lowdin_pairing_basis.zeroOverlapIndices(zero_spin)[0];
                    const auto co_density = lowdin_pairing_basis.coDensity(zero_overlap_index).component(zero_spin);

                    auto matrix_element = 0.0;

                    // Perform the contraction.
                    if (zero_spin == Spin::alpha) {
                        Tensor<Scalar, 0> element = operator_parameters_tensor_alpha.template einsum<2>("uv, vu ->", co_density.matrix());
                        matrix_element += element(0);
                    } else {
                        Tensor<Scalar, 0> element = operator_parameters_tensor_beta.template einsum<2>("uv, vu ->", co_density.matrix());
                        matrix_element += element(0);
                    }

                    evaluated_operator(i, j) += reduced_overlap * matrix_element;
                }
                // If there are two or more zero overlap values, the matrix element will be zero. No further if-clause is needed.
            }
        }

        // Return the matrix representation of the evaluated operator.
        return evaluated_operator;
    }


    /**
     * Evaluate any overlap operator in this non-orthogonal state basis.
     *
     * @return The matrix representation of this overlap operator expressed in the non-orthogonal state basis.
     *
     * @note This method requires no parameters, since the overlap in AO basis is a member of the non-orthogonal state basis.
     */
    Matrix evaluateOverlapOperator() const {

        // We start by initializing a zero matrix of the correct dimension.
        // Since the matrix is always squared, we only need one parameter to define its dimensions.
        Matrix evaluated_overlap = Matrix::Zero(this->numberOfBasisStates());

        // Loop over all basis states and calculate the correct overlap elements using the biorthogonalized basis of each combination of states.
        for (size_t i = 0; i < this->numberOfBasisStates(); i++) {
            for (size_t j = 0; j < this->numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N_a, this->N_b};

                // Fill in the overlaps in the new matrix representation in this non-orthogonal basis.
                evaluated_overlap(i, j) += lowdin_pairing_basis.totalOverlap();
            }
        }

        // Return the overlap matrix.
        return evaluated_overlap;
    }


    /**
     * Evaluate any scalar two-electron operator in this non-orthogonal state basis.
     *
     * @param g_op      The scalar two-electron operator to be evaluated.
     *
     * @return The matrix representation of this operator expressed in the non-orthogonal state basis.
     *
     * @note The given scalar two-electron operator must be expressed in the non-orthogonal AO basis.
     * @note This implementation is based on equation 65 from the 2021 paper by Hugh Burton (https://aip.scitation.org/doi/abs/10.1063/5.0045442).
     * @note This implementation only works for two electron operators where the mixed components are equal to zero.
     */
    Matrix evaluateTwoElectronOperator(const TwoElectronOperator& g_op) const {

        // We start by initializing a zero matrix of the correct dimension.
        // Since the matrix is always squared, we only need one parameter to define its dimensions.
        Matrix evaluated_operator = Matrix::Zero(this->numberOfBasisStates());

        // Loop over all basis states and calculate each matrix element using the generalized Slater-Condon rules.
        for (size_t i = 0; i < this->numberOfBasisStates(); i++) {
            for (size_t j = 0; j < this->numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N_a, this->N_b};

                // Check the number of zeros in the biorthogonal basis, as this influences how the matrix elements are calculated.
                const auto number_of_zeros = lowdin_pairing_basis.numberOfZeroOverlaps();

                // Use the correct formula, depending on how many zero overlaps there are.
                // If there are zero overlap values, we perform the following calculation.
                if (number_of_zeros == 0) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap and the weighted co-density matrix from the biorthogonal basis.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();
                    const auto weighted_co_density = lowdin_pairing_basis.weightedCoDensity();

                    // Perform the contractions. We need to perform two contractions (one for the direct component, one for the exchange component).
                    // In order to perserve readability of the contractions, and keep the exact link with the theory, we split each contraction in two.
                    Matrix intermediate_direct_contraction_1 = g_op.alphaAlpha().parameters().template einsum<2>("utvs, tu -> vs", weighted_co_density.alpha().matrix()).asMatrix();
                    Matrix intermediate_direct_contraction_2 = g_op.betaBeta().parameters().template einsum<2>("utvs, tu -> vs", weighted_co_density.beta().matrix()).asMatrix();

                    // The complete first direct contraction can the be written as follows.
                    Matrix direct_alpha_beta = intermediate_direct_contraction_1 + intermediate_direct_contraction_2;

                    // We return the representation to a tensor, in order to perform the final contractions.
                    Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> tensor_map {direct_alpha_beta.data(), direct_alpha_beta.rows(), direct_alpha_beta.cols()};
                    Tensor<Scalar, 2> direct_alpha_beta_tensor = Tensor<Scalar, 2>(tensor_map);

                    // We can now calculate the alpha and beta contributions with the next set of contractions.
                    Tensor<Scalar, 0> direct_element_a = direct_alpha_beta_tensor.template einsum<2>("vs, sv ->", weighted_co_density.alpha().matrix());
                    Tensor<Scalar, 0> direct_element_b = direct_alpha_beta_tensor.template einsum<2>("vs, sv ->", weighted_co_density.beta().matrix());

                    const auto direct_element = direct_element_a(0) + direct_element_b(0);

                    // Next, we calculate the exchange elements.
                    Tensor<Scalar, 2> intermediate_exchange_contraction_1 = g_op.alphaAlpha().parameters().template einsum<2>("utvs, su -> tv", weighted_co_density.alpha().matrix());
                    Tensor<Scalar, 0> exchange_element_a = intermediate_exchange_contraction_1.template einsum<2>("tv, tv ->", weighted_co_density.alpha().matrix());

                    Tensor<Scalar, 2> intermediate_exchange_contraction_2 = g_op.betaBeta().parameters().template einsum<2>("utvs, su -> tv", weighted_co_density.beta().matrix());
                    Tensor<Scalar, 0> exchange_element_b = intermediate_exchange_contraction_2.template einsum<2>("tv, tv ->", weighted_co_density.beta().matrix());

                    // We can now add the total contrinution to the corresponding metrix element.
                    evaluated_operator(i, j) += 0.5 * reduced_overlap * (direct_element - exchange_element_a(0) - exchange_element_b(0));
                }
                // If there is one zero overlap value, we perform the following calculation.
                else if (number_of_zeros == 1) {

                    // Determine whether the zero stems from the alpha or the beta component.
                    auto zero_spin = Spin::alpha;
                    auto non_zero_spin = Spin::beta;

                    if (lowdin_pairing_basis.numberOfZeroOverlaps(Spin::beta) > lowdin_pairing_basis.numberOfZeroOverlaps(Spin::alpha)) {
                        zero_spin = Spin::beta;
                        non_zero_spin = Spin::alpha;
                    }

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();

                    // Next, calculate the weighted co-density matrix and the co-density of the orbital corresponding to the zero overlap value.
                    // We know there is only one zero overlap orbital, so we acces the first index of the vector.
                    const auto zero_overlap_index = lowdin_pairing_basis.zeroOverlapIndices(zero_spin)[0];
                    const auto co_density = lowdin_pairing_basis.coDensity(zero_overlap_index).component(zero_spin);
                    const auto weighted_co_density = lowdin_pairing_basis.weightedCoDensity();

                    // Perform the contractions. We need to perform two contractions (one for the direct component, one for the exchange component).
                    // In order to perserve readability of the contractions, and keep the exact link with the theory, we split each contraction in two.
                    Tensor<Scalar, 2> intermediate_direct_contraction_1 = g_op.alphaAlpha().parameters().template einsum<2>("utvs, tu -> vs", co_density.matrix());
                    Tensor<Scalar, 0> direct_element_a = intermediate_direct_contraction_1.template einsum<2>("vs, sv ->", weighted_co_density.alpha().matrix());

                    Tensor<Scalar, 2> intermediate_direct_contraction_2 = g_op.betaBeta().parameters().template einsum<2>("utvs, tu -> vs", co_density.matrix());
                    Tensor<Scalar, 0> direct_element_b = intermediate_direct_contraction_2.template einsum<2>("vs, sv ->", weighted_co_density.beta().matrix());

                    auto direct_element = direct_element_a(0) + direct_element_b(0);

                    Tensor<Scalar, 2> intermediate_exchange_contraction = g_op.pureComponent(zero_spin).parameters().template einsum<2>("utvs, su -> tv", co_density.matrix());
                    Tensor<Scalar, 0> exchange_element = intermediate_exchange_contraction.template einsum<2>("tv, tv ->", weighted_co_density.component(zero_spin).matrix());

                    evaluated_operator(i, j) += (reduced_overlap * (direct_element - exchange_element(0)));
                }
                // If there are two zero overlap values, we perform the following calculation.
                else if (number_of_zeros == 2) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();

                    // First we merge the alpha and beta index vectors which denote the zero overlaps.
                    // Next, calculate the co-density of the orbitals corresponding to the zero overlap values.
                    const auto zero_overlap_index_1 = lowdin_pairing_basis.zeroOverlapIndices()[0].first;
                    const auto zero_overlap_index_2 = lowdin_pairing_basis.zeroOverlapIndices()[1].first;
                    const auto co_density_1 = lowdin_pairing_basis.coDensity(zero_overlap_index_1);

                    Matrix active_co_density = Matrix::Zero(co_density_1.alpha().matrix().dimension());

                    if (lowdin_pairing_basis.zeroOverlapIndices()[0].second == Spin::alpha) {
                        active_co_density += lowdin_pairing_basis.coDensity(zero_overlap_index_2).alpha().matrix();
                    } else {
                        active_co_density += lowdin_pairing_basis.coDensity(zero_overlap_index_2).beta().matrix();
                    }

                    // Perform the contractions. We need to perform two contractions (one for the direct component, one for the exchange component).
                    // In order to perserve readability of the contractions, and keep the exact link with the theory, we split each contraction in two.
                    Tensor<Scalar, 2> intermediate_direct_contraction_1 = g_op.alphaAlpha().parameters().template einsum<2>("utvs, tu -> vs", co_density_1.alpha().matrix());
                    Tensor<Scalar, 0> direct_element_a = intermediate_direct_contraction_1.template einsum<2>("vs, sv ->", active_co_density);

                    Tensor<Scalar, 2> intermediate_direct_contraction_2 = g_op.betaBeta().parameters().template einsum<2>("utvs, tu -> vs", co_density_1.beta().matrix());
                    Tensor<Scalar, 0> direct_element_b = intermediate_direct_contraction_2.template einsum<2>("vs, sv ->", active_co_density);

                    auto direct_element = direct_element_a(0) + direct_element_b(0);

                    auto exchange_element = 0.0;
                    if (lowdin_pairing_basis.zeroOverlapIndices()[0].second == Spin::alpha) {
                        Tensor<Scalar, 2> intermediate_exchange_contraction = g_op.pureComponent(Spin::alpha).parameters().template einsum<2>("utvs, su -> tv", co_density_1.alpha().matrix());
                        Tensor<Scalar, 0> element = intermediate_exchange_contraction.template einsum<2>("tv, tv ->", active_co_density);
                        exchange_element += element(0);
                    } else {
                        Tensor<Scalar, 2> intermediate_exchange_contraction = g_op.pureComponent(Spin::beta).parameters().template einsum<2>("utvs, su -> tv", co_density_1.beta().matrix());
                        Tensor<Scalar, 0> element = intermediate_exchange_contraction.template einsum<2>("tv, tv ->", active_co_density);
                        exchange_element += element(0);
                    }

                    evaluated_operator(i, j) += (reduced_overlap * (direct_element - exchange_element));
                }
                // If there are more than two zero overlap values, the matrix element will be zero. No further if-clause is needed.
            }
        }

        // Return the matrix representation of the evaluated operator.
        return evaluated_operator;
    }
};


/**
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _Scalar>
struct JacobiRotatableTraits<UNonOrthogonalStateBasis<_Scalar>> {

    // The type of Jacobi rotation that is naturally related to a `RNonOrthogonalStateBasis`.
    using JacobiRotationType = JacobiRotation;
};


/**
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _Scalar>
struct BasisTransformableTraits<UNonOrthogonalStateBasis<_Scalar>> {

    // The type of transformation that is naturally related to a `RNonOrthogonalStateBAsis`.
    using Transformation = UTransformation<_Scalar>;
};


}  // namespace GQCP
