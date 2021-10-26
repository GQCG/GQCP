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


#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Mathematical/Representation/Tensor.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Utilities/CRTP.hpp"
#include "Utilities/Eigen.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {

/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename NonOrthogonalStateBasis>
struct NonOrthogonalStateBasisTraits {};

/**
 *  A type that provides compile-time information on the operators associated with non-orthogonal state bases that is otherwise not accessible through a public class alias.
 */
template <typename NonOrthogonalStateBasis>
struct NOSBasisOperatorTraits {};

/*
 *  MARK: SimpleNonOrthogonalStateBasis
 */

/**
 *  A basis formed by any number of non-orthogonal states, in the form of `R/U/GTransformations.
 *
 *  @tparam _ExpansionScalar        The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
 */
template <typename _Scalar, typename _DerivedNonOrthogonalStateBasis>
class SimpleNonOrthogonalStateBasis:
    public CRTP<_DerivedNonOrthogonalStateBasis>,
    public BasisTransformable<_DerivedNonOrthogonalStateBasis>,
    public JacobiRotatable<_DerivedNonOrthogonalStateBasis> {

public:
    // The scalar type used to represent the expansion coefficients of the given non-orthogonal states: real or complex.
    using Scalar = _Scalar;

    // The type of the derived Lowdin pairing Basis from this parent class.
    using DerivedNonOrthogonalStateBasis = _DerivedNonOrthogonalStateBasis;

    // The type of matrix associated with this kind of NonOrthogonalStateBasis.
    using Matrix = SquareMatrix<Scalar>;

    // The type of transformation that is naturally related to a `NonOrthogonalStateBasis`. Can be R/U/GTransformation.
    using Transformation = typename NonOrthogonalStateBasisTraits<DerivedNonOrthogonalStateBasis>::Transformation;

    // The biorthogonal basis related to the basis states in this type of non-orthogonal basis (R/U/GLowdinPairingBasis).
    using BiorthogonalBasis = typename NonOrthogonalStateBasisTraits<DerivedNonOrthogonalStateBasis>::BiorthogonalBasis;

    // The type of Jacobi rotation that is naturally related to the derived non orthogonal-basis.
    using JacobiRotationType = typename JacobiRotatableTraits<DerivedNonOrthogonalStateBasis>::JacobiRotationType;

    // The second-quantized representation of the Hamiltonian that can be evaluated in this basis.
    using Hamiltonian = typename NOSBasisOperatorTraits<DerivedNonOrthogonalStateBasis>::Hamiltonian;

    // The scalar one-electron operator that can be evaluated by the non-orthogonal basis.
    using OneElectronOperator = typename NOSBasisOperatorTraits<DerivedNonOrthogonalStateBasis>::OneElectronOperator;

    // The scalar two-electron operator that can be evaluated by the non-orthogonal basis.
    using TwoElectronOperator = typename NOSBasisOperatorTraits<DerivedNonOrthogonalStateBasis>::TwoElectronOperator;

    // The vector containing the basis state of the associated type of transformations.
    using States = std::vector<Transformation>;


protected:
    // The vector containing the non-orthogonal basis states.
    States basis_states;

    // The overlap operator in AO basis, constructed from the spinor/spin-orbital basis.
    OneElectronOperator overlap_operator_AO;

    // The total number of occupied orbitals.
    size_t N;

    // The threshold used to determine zero values.
    double zero_threshold;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `SimpleNonOrthogonalStateBasis` from any number of non orthogonal states.
     *
     *  @param basis_state_vector               The vector containing the non-orthogonal basis states.
     *  @param S_AO                             The overlap operator in AO basis, constructed from the spinor/spin-orbital used to calculate the non-orthogonal states.
     *  @param number_of_occupied_orbitals      The total number of occupied orbitals in the system.
     *  @param threshold                        The threshold at which a value is verified to be zero or not. The default is 1e-8.
     */
    SimpleNonOrthogonalStateBasis<Scalar, DerivedNonOrthogonalStateBasis>(const States& basis_state_vector, const OneElectronOperator& S_AO, const size_t number_of_occupied_orbitals, const double threshold = 1e-8) :
        basis_states {basis_state_vector},
        overlap_operator_AO {S_AO},
        N {number_of_occupied_orbitals},
        zero_threshold {threshold} {

        // The basis states must have the same dimensions.
        for (size_t i = 0; i < basis_state_vector.size(); i++) {
            if (basis_state_vector[0].dimension() != basis_state_vector[i].dimension()) {
                throw std::invalid_argument("NonOrthogonalStateBasis<Scalar>(const States& basis_state_vector, const OneElectronOperator& S_AO, const size_t number_of_occupied_orbitals, const double threshold = 1e-8): The given basis states do not have the same dimensions.");
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
     * @note We return the dimension of the first state, as the constructor checks that all states have the same dimension.
     */
    const size_t basisStateDimension() const { return this->basis_states[0].dimension(); }

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


    /*
     *  MARK: Conforming to `BasisTransformable`
     */

    /**
     *  Apply the basis transformation and return the result.
     *
     *  @param T            The basis transformation.
     *
     *  @return The basis-transformed object.
     */
    DerivedNonOrthogonalStateBasis transformed(const Transformation& T) const override {

        auto result = this->derived();

        for (size_t i = 0; i < result.numberOfBasisStates(); i++) {
            result.basis_states[i].transform(T);
        }

        return result;
    }

    // Allow the `rotate` method from `BasisTransformable`, since there's also a `rotate` from `JacobiRotatable`.
    using BasisTransformable<DerivedNonOrthogonalStateBasis>::rotate;

    // Allow the `rotated` method from `BasisTransformable`, since there's also a `rotated` from `JacobiRotatable`.
    using BasisTransformable<DerivedNonOrthogonalStateBasis>::rotated;


    /*
     *  MARK: Conforming to `JacobiRotatable`.
     */

    /**
     *  Apply the Jacobi rotation and return the result.
     *
     *  @param jacobi_rotation          The Jacobi rotation.
     *
     *  @return The Jacobi-rotated object.
     */
    DerivedNonOrthogonalStateBasis rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto J = Transformation::FromJacobi(jacobi_rotation, this->basisStateDimension());
        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<DerivedNonOrthogonalStateBasis>::rotate;


    /**
     * MARK: Operator evaluations
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
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> operator_tensor_map {f_op.parameters().data(), f_op.parameters().rows(), f_op.parameters().cols()};
        Tensor<Scalar, 2> operator_parameters_tensor = Tensor<Scalar, 2>(operator_tensor_map);

        // Loop over all basis states and calculate each matrix element using the generalized Slater-Condon rules.
        for (size_t i = 0; i < this->numberOfBasisStates(); i++) {
            for (size_t j = 0; j < this->numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N};

                // Check the number of zeros in the biorthogonal basis, as this influences how the matrix elements are calculated.
                const auto number_of_zeros = lowdin_pairing_basis.numberOfZeroOverlaps();

                // Use the correct formula, depending on how many zero overlaps there are.
                // If there are zero overlap values, we perform the following calculation.
                if (number_of_zeros == 0) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap and the weighted co-density matrix from the biorthogonal basis.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();
                    const auto weighted_co_density = lowdin_pairing_basis.weightedCoDensity();

                    // Perform the contraction.
                    Tensor<Scalar, 0> matrix_element = reduced_overlap * operator_parameters_tensor.template einsum<2>("uv, vu ->", weighted_co_density.matrix());
                    evaluated_operator(i, j) += matrix_element(0);
                }
                // If there is one zero overlap value, we perform the following calculation.
                else if (number_of_zeros == 1) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();

                    // Next, calculate the co-density matrix of the orbital corresponding to the zero overlap value.
                    // We know there is only one zero overlap orbital, so we acces the first index of the vector.
                    const auto zero_overlap_index = lowdin_pairing_basis.zeroOverlapIndices()[0];
                    const auto co_density = lowdin_pairing_basis.coDensity(zero_overlap_index);

                    // Perform the contraction.
                    Tensor<Scalar, 0> matrix_element = reduced_overlap * operator_parameters_tensor.template einsum<2>("uv, vu ->", co_density.matrix());
                    evaluated_operator(i, j) += matrix_element(0);
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
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N};

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
     */
    Matrix evaluateTwoElectronOperator(const TwoElectronOperator& g_op) const {

        // We start by initializing a zero matrix of the correct dimension.
        // Since the matrix is always squared, we only need one parameter to define its dimensions.
        Matrix evaluated_operator = Matrix::Zero(this->numberOfBasisStates());

        // Loop over all basis states and calculate each matrix element using the generalized Slater-Condon rules.
        for (size_t i = 0; i < this->numberOfBasisStates(); i++) {
            for (size_t j = 0; j < this->numberOfBasisStates(); j++) {

                // The first step is to create a biorthogonal basis from the two states that are being looped over.
                const BiorthogonalBasis lowdin_pairing_basis {this->basisState(i), this->basisState(j), this->overlap_operator_AO, this->N};

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
                    Tensor<Scalar, 2> intermediate_contraction_1 = g_op.parameters().template einsum<2>("utvs, tu -> vs", weighted_co_density.matrix());
                    Tensor<Scalar, 0> direct_element = intermediate_contraction_1.template einsum<2>("vs, sv ->", weighted_co_density.matrix());

                    Tensor<Scalar, 2> intermediate_contraction_2 = g_op.parameters().template einsum<2>("utvs, su -> tv", weighted_co_density.matrix());
                    Tensor<Scalar, 0> exchange_element = intermediate_contraction_2.template einsum<2>("tv, tv ->", weighted_co_density.matrix());

                    evaluated_operator(i, j) += (0.5 * reduced_overlap * (direct_element(0) - exchange_element(0)));
                }
                // If there is one zero overlap value, we perform the following calculation.
                else if (number_of_zeros == 1) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();

                    // Next, calculate the weighted co-density matrix and the co-density of the orbital corresponding to the zero overlap value.
                    // We know there is only one zero overlap orbital, so we acces the first index of the vector.
                    const auto zero_overlap_index = lowdin_pairing_basis.zeroOverlapIndices()[0];
                    const auto co_density = lowdin_pairing_basis.coDensity(zero_overlap_index);
                    const auto weighted_co_density = lowdin_pairing_basis.weightedCoDensity();

                    // Perform the contractions. We need to perform two contractions (one for the direct component, one for the exchange component).
                    // In order to perserve readability of the contractions, and keep the exact link with the theory, we split each contraction in two.
                    Tensor<Scalar, 2> intermediate_contraction_1 = g_op.parameters().template einsum<2>("utvs, tu -> vs", co_density.matrix());
                    Tensor<Scalar, 0> direct_element = intermediate_contraction_1.template einsum<2>("vs, sv ->", weighted_co_density.matrix());

                    Tensor<Scalar, 2> intermediate_contraction_2 = g_op.parameters().template einsum<2>("utvs, su -> tv", co_density.matrix());
                    Tensor<Scalar, 0> exchange_element = intermediate_contraction_2.template einsum<2>("tv, tv ->", weighted_co_density.matrix());

                    evaluated_operator(i, j) += (reduced_overlap * (direct_element(0) - exchange_element(0)));
                }
                // If there are two zero overlap values, we perform the following calculation.
                else if (number_of_zeros == 2) {

                    // In order to calcluate the matrix element when there are no zero overlap values, we need the reduced overlap.
                    const auto reduced_overlap = lowdin_pairing_basis.reducedOverlap();

                    // Next, calculate the co-density of the orbitals corresponding to the zero overlap values.
                    // We know there are only two zero overlap orbitals, so we acces the first and second index of the vector.
                    const auto zero_overlap_index_1 = lowdin_pairing_basis.zeroOverlapIndices()[0];
                    const auto zero_overlap_index_2 = lowdin_pairing_basis.zeroOverlapIndices()[1];
                    const auto co_density_1 = lowdin_pairing_basis.coDensity(zero_overlap_index_1);
                    const auto co_density_2 = lowdin_pairing_basis.coDensity(zero_overlap_index_2);

                    // Perform the contractions. We need to perform two contractions (one for the direct component, one for the exchange component).
                    // In order to perserve readability of the contractions, and keep the exact link with the theory, we split each contraction in two.
                    Tensor<Scalar, 2> intermediate_contraction_1 = g_op.parameters().template einsum<2>("utvs, tu -> vs", co_density_1.matrix());
                    Tensor<Scalar, 0> direct_element = intermediate_contraction_1.template einsum<2>("vs, sv ->", co_density_2.matrix());

                    Tensor<Scalar, 2> intermediate_contraction_2 = g_op.parameters().template einsum<2>("utvs, su -> tv", co_density_1.matrix());
                    Tensor<Scalar, 0> exchange_element = intermediate_contraction_2.template einsum<2>("tv, tv ->", co_density_2.matrix());

                    evaluated_operator(i, j) += (reduced_overlap * (direct_element(0) - exchange_element(0)));
                }
                // If there are more than two zero overlap values, the matrix element will be zero. No further if-clause is needed.
            }
        }

        // Return the matrix representation of the evaluated operator.
        return evaluated_operator;
    }
};


}  // namespace GQCP
