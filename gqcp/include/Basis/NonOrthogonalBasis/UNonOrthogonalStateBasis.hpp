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
#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
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
    OneElectronOperator overlap_operator_AO;

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
    UNonOrthogonalStateBasis<Scalar>(const States& basis_state_vector, const OneElectronOperator& S_AO, const size_t number_of_occupied_alpha_orbitals, const size_t number_of_occupied_beta_orbitals, const double threshold = 1e-8) :
        basis_states {basis_state_vector},
        overlap_operator_AO {S_AO},
        N_a {number_of_occupied_alpha_orbitals},
        N_b {number_of_occupied_beta_orbitals},
        zero_threshold {threshold} {

        // The basis states must have the same dimensions.
        for (size_t i = 0; i < basis_state_vector.size(); i++) {
            if (basis_state_vector[0].alpha().dimension() != basis_state_vector[i].alpha().dimension()) {
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
    Self rotated(const JacobiRotationType& jacobi_rotation) const override {

        const auto Ja = UTransformationComponent<Scalar>::FromJacobi(jacobi_rotation, this->basisStateDimension());
        const auto Jb = UTransformationComponent<Scalar>::FromJacobi(jacobi_rotation, this->basisStateDimension());
        Transformation J {Ja, Jb};

        return this->rotated(J);
    }

    // Allow the `rotate` method from `JacobiRotatable`, since there's also a `rotate` from `BasisTransformable`.
    using JacobiRotatable<Self>::rotate;
};


/*
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


/*
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
