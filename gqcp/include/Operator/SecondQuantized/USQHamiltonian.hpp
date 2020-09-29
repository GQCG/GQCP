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


#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QuantumChemical/Spin.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A class for representing the electronic Hamiltonian, expressed in an unrestricted spinor basis.
 * 
 *  @tparam Scalar      the scalar type of the second-quantized parameters (i.e. the integrals)
 */
template <typename Scalar>
class USQHamiltonian {
private:
    std::array<SQHamiltonian<Scalar>, 2> sq_hamiltonians;  // array that holds the individual SQHamiltonian for the pure alpha and beta components (in that order)

    std::vector<ScalarSQTwoElectronOperator<Scalar>> two_op_mixed;  // the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
    ScalarSQTwoElectronOperator<Scalar> total_two_op_mixed;         // the total alpha & beta mixed two-electron operator (whose integrals are represented as g_aabb)

public:
    /*
     *  CONSTRUCTORS
     */
    USQHamiltonian() = default;

    /**
     *  @param sq_hamiltonian_alpha      the alpha Hamiltonian
     *  @param sq_hamiltonian_beta       the beta Hamiltonian
     *  @param two_op_mixed              the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
     */
    USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const std::vector<ScalarSQTwoElectronOperator<Scalar>>& two_op_mixed) :
        sq_hamiltonians {sq_hamiltonian_alpha, sq_hamiltonian_beta},
        two_op_mixed {two_op_mixed} {

        // Check if the dimensions are compatible
        const auto dim = sq_hamiltonians[Spin::alpha].numberOfOrbitals();

        if (sq_hamiltonians[Spin::beta].numberOfOrbitals() != dim) {
            throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed): The dimensions of the alpha and beta Hamiltonian are incompatible");
        }

        for (const auto& two_op : this->two_op_mixed) {
            if (two_op.numberOfOrbitals() != dim) {
                throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed): The dimensions of the mixed two electron operator are incompatible with the Hamiltonian");
            }
        }

        // Calculate the total two-electron operator
        QCRankFourTensor<Scalar> total_two_op_par {dim};
        total_two_op_par.setZero();
        for (const auto& two_op : this->two_op_mixed) {
            total_two_op_par += two_op.parameters().Eigen();
        }
        this->total_two_op_mixed = ScalarSQTwoElectronOperator<Scalar>(total_two_op_par);
    }

    /**
     *  @param sq_hamiltonian_alpha      the alpha Hamiltonian
     *  @param sq_hamiltonian_beta       the beta Hamiltonian
     *  @param two_op_mixed              the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
     */
    USQHamiltonian(const SQHamiltonian<Scalar>& sq_hamiltonian_alpha, const SQHamiltonian<Scalar>& sq_hamiltonian_beta, const ScalarSQTwoElectronOperator<Scalar>& two_op_mixed) :
        USQHamiltonian(sq_hamiltonian_alpha, sq_hamiltonian_beta,
                       std::vector<ScalarSQTwoElectronOperator<Scalar>> {two_op_mixed}) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Construct the molecular Hamiltonian in a given single-particle basis
     *
     *  @param spinor_basis         the initial unrestricted spinor basis in which the components of the Hamiltonian should be expressed
     *  @param molecule             the molecule on which the  spinor basis is based
     *
     *  @return a second-quantized molecular unrestricted Hamiltonian. The molecular unrestricted Hamiltonian has:
     *      - restricted Hamiltonians for the alpha and beta components
     *      - mixed two-electron contributions
     *
     *  Note that this named constructor is only available for real matrix representations
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> Molecular(const USpinorBasis<Z, GTOShell>& u_spinor_basis, const Molecule& molecule) {

        const SQHamiltonian<Scalar> sq_hamiltonian_alpha = SQHamiltonian<double>::Molecular(u_spinor_basis.spinorBasis(Spin::alpha), molecule);
        const SQHamiltonian<Scalar> sq_hamiltonian_beta = SQHamiltonian<double>::Molecular(u_spinor_basis.spinorBasis(Spin::beta), molecule);

        // Initial basis for alpha and beta are identical so the mixed integrals are identical to spin specific components
        const ScalarSQTwoElectronOperator<double> two_op_mixed = sq_hamiltonian_alpha.twoElectron();

        return USQHamiltonian(sq_hamiltonian_alpha, sq_hamiltonian_beta, two_op_mixed);
    }

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return if the alpha and beta components of the unrestricted Hamiltonian are of the same dimension
     */
    bool areSpinHamiltoniansOfSameDimension() const { return this->spinHamiltonian(Spin::alpha).numberOfOrbitals() == this->spinHamiltonian(Spin::beta).numberOfOrbitals(); }

    /**
     *  Constrain a spin component of the unrestricted Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_electron_op           the one-electron operator used as a constraint
     *  @param lambda                    the Lagrangian multiplier for the constraint
     *  @param sigma                     the spin sigma component
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> constrain(const ScalarSQOneElectronOperator<double>& one_electron_op, const double lambda, const Spin& sigma) const {

        auto const constrained_component = this->sq_hamiltonians[sigma] - lambda * one_electron_op;

        if (sigma == Spin::beta) {
            return USQHamiltonian(this->sq_hamiltonians[Spin::alpha], constrained_component, this->two_op_mixed);
        } else {
            return USQHamiltonian(constrained_component, this->sq_hamiltonians[Spin::beta], this->two_op_mixed);
        }
    }


    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t numberOfOrbitals() const { return this->sq_hamiltonians[Spin::alpha].numberOfOrbitals() + this->sq_hamiltonians[Spin::beta].numberOfOrbitals(); }

    /**
     *  In-place rotate the matrix representations of the Hamiltonian
     *      
     *  @param U    the unitary rotation matrix between the old and the new orbital basis
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        this->sq_hamiltonians[Spin::alpha].rotate(U);
        this->sq_hamiltonians[Spin::beta].rotate(U);
    }


    /**
     *  In-place rotate the matrix representation of one of the spin components of the Hamiltonian
     *
     *  @param U                    the unitary rotation matrix between the old and the new orbital basis
     *  @param sigma                the spin sigma component
     */
    void rotate(const TransformationMatrix<Scalar>& U, const Spin& sigma) { this->transform(U, sigma); }


    /**
     *  In-place rotate the matrix representations of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        this->sq_hamiltonians[Spin::alpha].rotate(jacobi_rotation_parameters);
        this->sq_hamiltonians[Spin::beta].rotate(jacobi_rotation_parameters);
    }


    /**
     *  In-place rotate the matrix representation of one of the spin components of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     *  @param sigma                            the spin sigma component
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters, const Spin& sigma) {
        this->transform(jacobi_rotation_parameters, sigma);
    }


    /**
     *  @param sigma                    the spin sigma component
     * 
     *  @return the pure contributions of the requested component of the unrestricted Hamiltonian 
     */
    const SQHamiltonian<Scalar>& spinHamiltonian(const Spin& sigma) const { return this->sq_hamiltonians[sigma]; }

    /**
     *  In-place transform the matrix representations of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        this->sq_hamiltonians[Spin::alpha].transform(T);
        this->sq_hamiltonians[Spin::beta].transform(T);

        // Transform the mixed
        this->total_two_op_mixed.transform(T);
        for (auto& two_op : this->two_op_mixed) {
            two_op.transform(T);
        }
    }


    /**
     *  In-place transform the matrix representations of a single spin component of the unrestricted Hamiltonian
     * 
     *  @param T                        the transformation matrix between the old and the new orbital basis
     *  @param component                the spin component
     */
    void transform(const TransformationMatrix<Scalar>& T, const Spin& component) {

        this->sq_hamiltonians[component].transform(T);

        // Transform the mixed two-electron operators and their total
        auto new_two_electron_parameters = this->total_two_op_mixed.parameters();

        // transform the two electron parameters "g_aabb" to "g_a'a'bb" when alpha is chosen or to "g_aab'b'" for beta.
        const size_t first_contraction_index = 2 * component;
        const size_t second_contraction_index = 2 * component + 1;
        new_two_electron_parameters.template contractWithMatrix<Scalar>(T, first_contraction_index);
        new_two_electron_parameters.template contractWithMatrix<Scalar>(T, second_contraction_index);
        this->total_two_op_mixed = ScalarSQTwoElectronOperator<Scalar>(new_two_electron_parameters);

        for (auto& two_op : this->two_op_mixed) {

            auto new_two_electron_parameters = two_op.parameters();
            // transform the two electron parameters "g_aabb" to "g_a'a'bb"
            new_two_electron_parameters.template contractWithMatrix<Scalar>(T, first_contraction_index);
            new_two_electron_parameters.template contractWithMatrix<Scalar>(T, second_contraction_index);
            two_op = ScalarSQTwoElectronOperator<Scalar>(new_two_electron_parameters);
        }
    }


    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const std::vector<ScalarSQTwoElectronOperator<Scalar>>& twoElectronContributionsMixed() const { return this->two_op_mixed; }

    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectronMixed() const { return this->total_two_op_mixed; }
};


}  // namespace GQCP
