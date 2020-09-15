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


#include "Basis/SpinorBasis/Spin.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
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
    ScalarUSQOneElectronOperator<Scalar> total_one_op;  // one-electron interactions (i.e. the core Hamiltonian)
    ScalarUSQTwoElectronOperator<Scalar> total_two_op;  // two-electron interactions

    std::vector<ScalarUSQOneElectronOperator<Scalar>> one_ops;  // the core (i.e. one-electron) contributions to the Hamiltonian
    std::vector<ScalarUSQTwoElectronOperator<Scalar>> two_ops;  // the two-electron contributions to the Hamiltonian

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
    USQHamiltonian(const std::vector<ScalarUSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarUSQTwoElectronOperator<Scalar>>& two_ops) :
        one_ops {one_ops},
        two_ops {two_ops} {

        // Check if the dimensions are compatible
        const std::invalid_argument dimension_error("USQHamiltonian::USQHamiltonian(const std::vector<ScalarUSQOneElectronOperator<Scalar>& one_ops, const std::vector<ScalarUSQTwoElectronOperator<Scalar>& two_ops: The dimensions of the operators and coefficients matrix are incompatible");

        // Start with the alpha and beta dimensions of the one-electron parts
        const auto dim_a = one_ops[0].parameters(Spin::alpha).numberOfOrbitals();
        const auto dim_b = one_ops[0].parameters(Spin::beta).numberOfOrbitals();

        for (const auto& one_op : this->one_ops) {
            if ((one_op.parameters(Spin::alpha).numberOfOrbitals() != dim_a) || (one_op.parameters(Spin::beta).numberOfOrbitals() != dim_b)) {
                throw dimension_error;
            }
        }

        // Do the same for all four two electron components
        const auto dim_aa = two_ops[0].parameters(Spin::alpha, Spin::alpha).numberOfOrbitals();
        const auto dim_ab = two_ops[0].parameters(Spin::alpha, Spin::beta).numberOfOrbitals();
        const auto dim_ba = two_ops[0].parameters(Spin::beta, Spin::alpha).numberOfOrbitals();
        const auto dim_bb = two_ops[0].parameters(Spin::beta, Spin::beta).numberOfOrbitals();

        for (const auto& two_op : this->two_ops) {
            if ((two_op.parameters.parameters(Spin::alpha, Spin::alpha).numberOfOrbitals() != dim_aa) || (two_op.parameters.parameters(Spin::alpha, Spin::beta).numberOfOrbitals() != dim_ab) || (two_op.parameters.parameters(Spin::beta, Spin::alpha).numberOfOrbitals() != dim_ba) || (two_op.parameters.parameters(Spin::beta, Spin::beta).numberOfOrbitals() != dim_bb)) {
                throw dimension_error;
            }
        }

        // Calculate the total one-electron operator
        QCMatrix<Scalar> total_one_op_par_a {dim_a};
        QCMatrix<Scalar> total_one_op_par_b {dim_b};

        total_one_op_par_a.setZero();
        total_one_op_par_b.setZero();

        for (const auto& one_op : this->one_ops) {
            total_one_op_par_a += one_op.parameters(Spin::alpha);
            total_one_op_par_b += one_op.parameters(Spin::beta);
        }

        this->total_one_op = ScalarUSQOneElectronOperator<Scalar>(total_one_op_par_a, total_one_op_par_b);

        // Calculate the total two-electron operator
        QCRankFourTensor<Scalar> total_two_op_par_aa(dim_aa);
        QCRankFourTensor<Scalar> total_two_op_par_ab(dim_ab);
        QCRankFourTensor<Scalar> total_two_op_par_ba(dim_ba);
        QCRankFourTensor<Scalar> total_two_op_par_bb(dim_bb);

        total_two_op_par_aa.setZero();
        total_two_op_par_ab.setZero();
        total_two_op_par_ba.setZero();
        total_two_op_par_bb.setZero();

        for (const auto& two_op : this->two_ops) {
            total_two_op_par_aa += two_op.parameters(Spin::alpha, Spin::alpha).Eigen();
            total_two_op_par_ab += two_op.parameters(Spin::alpha, Spin::beta).Eigen();
            total_two_op_par_ba += two_op.parameters(Spin::beta, Spin::alpha).Eigen();
            total_two_op_par_bb += two_op.parameters(Spin::beta, Spin::beta).Eigen();
        }

        this->total_two_op = ScalarUSQTwoElectronOperator<Scalar>(total_two_op_par_aa, total_two_op_par_ab, total_two_op_par_ba, total_two_op_par_bb);
    }


    /**
     *  @param sq_hamiltonian_alpha      the alpha Hamiltonian
     *  @param sq_hamiltonian_beta       the beta Hamiltonian
     *  @param two_op_mixed              the alpha & beta mixed two-electron operators (whose integrals are represented as g_aabb)
     */
    USQHamiltonian(const ScalarUSQOneElectronOperator<Scalar>& h, const ScalarUSQTwoElectronOperator<Scalar>& g) :
        USQHamiltonian(std::vector<ScalarUSQOneElectronOperator<Scalar>> {h}, std::vector<ScalarUSQTwoElectronOperator<Scalar>> {g}) {}


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

        // Create the alpha and beta one-electron operators separately and then encapsulate them
        const auto T_a = u_spinor_basis.quantize(Operator::Kinetic(), Spin::alpha);
        const auto V_a = u_spinor_basis.quantize(Operator::NuclearAttraction(molecule), Spin::alpha);
        ScalarSQOneElectronOperator<double> H_a = T_a + V_a;

        const auto T_b = u_spinor_basis.quantize(Operator::Kinetic(), Spin::beta);
        const auto V_b = u_spinor_basis.quantize(Operator::NuclearAttraction(molecule), Spin::beta);
        ScalarSQOneElectronOperator<double> H_b = T_b + V_b;

        const ScalarUSQOneElectronOperator<double> H {H_a.parameters(), H_b.parameters()};

        // Do the same for the two-electron operators
        const auto g_aa = u_spinor_basis.quantize(Operator::Coulomb(), Spin::alpha);
        const auto g_bb = u_spinor_basis.quantize(Operator::Coulomb(), Spin::beta);

        // Initial basis for alpha and beta are identical so the mixed integrals are identical to spin specific components
        // It would be better to allow the USpinorBasis to quantize these mixed components, as I'm not sure this method is correct.
        const auto g_ab = g_aa;
        const auto g_ba = g_bb;

        const ScalarUSQTwoElectronOperator<double> g {g_aa.parameters(), g_ab.parameters(), g_ba.parameters(), g_bb.parameters()};

        return USQHamiltonian(H, g);
    }

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return if the alpha and beta components of the unrestricted Hamiltonian are of the same dimension
     */
    bool areSpinHamiltoniansOfSameDimension() const { return this->spinHamiltonian(Spin::alpha).numberOfOrbitals() == this->spinHamiltonian(Spin::beta).numberOfOrbitals(); }

    // /**
    //  *  Constrain a spin component of the unrestricted Hamiltonian according to the convention: - lambda * constraint
    //  *
    //  *  @param one_electron_op           the one-electron operator used as a constraint
    //  *  @param lambda                    the Lagrangian multiplier for the constraint
    //  *  @param sigma                     the spin sigma component
    //  *
    //  *  @return a copy of the constrained Hamiltonian
    //  *
    //  *  Note that this method is only available for real matrix representations
    //  */
    // template <typename Z = Scalar>
    // enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> constrain(const ScalarSQOneElectronOperator<double>& one_electron_op, const double lambda, const Spin& sigma) const {

    //     auto const constrained_component = this->sq_hamiltonians[sigma] - lambda * one_electron_op;

    //     if (sigma == Spin::beta) {
    //         return USQHamiltonian(this->sq_hamiltonians[Spin::alpha], constrained_component, this->two_op_mixed);
    //     } else {
    //         return USQHamiltonian(constrained_component, this->sq_hamiltonians[Spin::beta], this->two_op_mixed);
    //     }
    // }


    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t numberOfOrbitals() const { return this->spinHamiltonian(Spin::alpha).numberOfOrbitals() + this->spinHamiltonian(Spin::beta).numberOfOrbitals(); }

    /**
     *  In-place rotate the matrix representations of the Hamiltonian
     *      
     *  @param U    the unitary rotation matrix between the old and the new orbital basis
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Rotate the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.rotate(U);
        }

        // Rotate the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.rotate(U);
        }

        // Rotate the totals
        this->total_one_op.rotate(U);
        this->total_two_op.rotate(U);
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

        // Transform the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.rotate(jacobi_rotation_parameters);
        }

        // Transform the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.rotate(jacobi_rotation_parameters);
        }

        // Transform the totals
        this->total_one_op.rotate(jacobi_rotation_parameters);
        this->total_two_op.rotate(jacobi_rotation_parameters);
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
     *  @return the pure contributions of the requested component of the unrestricted Hamiltonian *NOT YET CORRECT, DOESN'T CONTAIN VECTORS OF CONTRIBUTIONS*
     */
    const SQHamiltonian<Scalar> spinHamiltonian(const Spin& sigma) const {

        std::vector<ScalarSQOneElectronOperator<Scalar>> one_ops_sigma;
        std::vector<ScalarSQTwoElectronOperator<Scalar>> two_ops_sigma;

        switch (sigma) {
        case Spin::alpha: {
            for (auto& one_op : this->one_ops) {
                one_ops_sigma.push_back(ScalarSQOneElectronOperator<Scalar>(one_op.parameters(Spin::alpha)));
            }
            for (auto& two_op : this->two_ops) {
                two_ops_sigma.push_back(ScalarSQTwoElectronOperator<Scalar>(two_op.parameters(Spin::alpha, Spin::alpha)));
            }

            return SQHamiltonian<Scalar> {one_ops_sigma, two_ops_sigma};
            break;
        }

        case Spin::beta: {
            for (auto& one_op : this->one_ops) {
                one_ops_sigma.push_back(ScalarSQOneElectronOperator<Scalar>(one_op.parameters(Spin::beta)));
            }
            for (auto& two_op : this->two_ops) {
                two_ops_sigma.push_back(ScalarSQTwoElectronOperator<Scalar>(two_op.parameters(Spin::beta, Spin::beta)));
            }

            return SQHamiltonian<Scalar> {one_ops_sigma, two_ops_sigma};
            break;
        }
        }
    }

    /**
     *  In-place transform the matrix representations of the unrestricted Hamiltonian
     *
     *  @param T    the transformation matrix between the old and the new orbital basis
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.transform(T);
        }

        // Transform the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.transform(T);
        }

        // Transform the totals
        this->total_one_op.transform(T);
        this->total_two_op.transform(T);
    }


    /**
     *  In-place transform the matrix representations of a single spin component of the unrestricted Hamiltonian
     * 
     *  @param T                        the transformation matrix between the old and the new orbital basis
     *  @param sigma                    the spin sigma component
     */
    void transform(const TransformationMatrix<Scalar>& T, const Spin& sigma) {

        // initialize new total one-electron operators and set them to zero
        QCMatrix<Scalar> total_one_op_par_a {this->one_ops[0].parameters(Spin::alpha).numberOfOrbitals()};
        QCMatrix<Scalar> total_one_op_par_b {this->one_ops[0].parameters(Spin::beta).numberOfOrbitals()};
        total_one_op_par_a.setZero();
        total_one_op_par_b.setZero();

        // initialize new total two-electron operators and set them to zero
        QCRankFourTensor<Scalar> total_two_op_par_aa(two_ops[0].parameters(Spin::alpha, Spin::alpha).numberOfOrbitals());
        QCRankFourTensor<Scalar> total_two_op_par_ab(two_ops[0].parameters(Spin::alpha, Spin::beta).numberOfOrbitals());
        QCRankFourTensor<Scalar> total_two_op_par_ba(two_ops[0].parameters(Spin::beta, Spin::alpha).numberOfOrbitals());
        QCRankFourTensor<Scalar> total_two_op_par_bb(two_ops[0].parameters(Spin::beta, Spin::beta).numberOfOrbitals());

        total_two_op_par_aa.setZero();
        total_two_op_par_ab.setZero();
        total_two_op_par_ba.setZero();
        total_two_op_par_bb.setZero();

        switch (sigma) {
        case Spin::alpha: {

            // transform the alpha components of the one electron operators
            for (auto& one_op : this->one_ops) {
                one_op.parameters(Spin::alpha).basisTransform(T);
                total_one_op_par_a += one_op.parameters(Spin::alpha);
                total_one_op_par_b += one_op.parameters(Spin::beta);
            }
            this->total_one_op = ScalarUSQOneElectronOperator<Scalar>(total_one_op_par_a, total_one_op_par_b);

            // Transform the pure alpha and alpha-beta components of the two-electron operators when alpha is chosen
            // transform the two electron parameters "g_aabb" to "g_a'a'bb" when alpha is chosen
            for (auto& two_op : this->two_ops) {
                two_op.parameters(Spin::alpha, Spin::alpha).bassisTransform(T);
                two_op.parameters(Spin::alpha, Spin::beta).basisTransform(T);
                total_two_op_par_aa += two_op.parameters(Spin::alpha, Spin::alpha).Eigen();
                total_two_op_par_ab += two_op.parameters(Spin::alpha, Spin::beta).Eigen();
                total_two_op_par_ba += two_op.parameters(Spin::beta, Spin::alpha).Eigen();
                total_two_op_par_bb += two_op.parameters(Spin::beta, Spin::beta).Eigen();
            }
            this->total_two_op = ScalarUSQTwoElectronOperator<Scalar>(total_two_op_par_aa, total_two_op_par_ab, total_two_op_par_ba, total_two_op_par_bb);
            break;
        }

        case Spin::beta: {

            // transform the beta components of the one electron operators
            for (auto& one_op : this->one_ops) {
                one_op.parameters(Spin::beta).basisTransform(T);
                total_one_op_par_a += one_op.parameters(Spin::alpha);
                total_one_op_par_b += one_op.parameters(Spin::beta);
            }
            this->total_one_op = ScalarUSQOneElectronOperator<Scalar>(total_one_op_par_a, total_one_op_par_b);

            // Transform the pure beta and beta-alpha components of the two-electron operators when beta is chosen
            // transform the two electron parameters "g_aab'b'" for beta.
            for (auto& two_op : this->two_ops) {
                two_op.parameters(Spin::beta, Spin::beta).bassisTransform(T);
                two_op.parameters(Spin::beta, Spin::alpha).basisTransform(T);
                total_two_op_par_aa += two_op.parameters(Spin::alpha, Spin::alpha).Eigen();
                total_two_op_par_ab += two_op.parameters(Spin::alpha, Spin::beta).Eigen();
                total_two_op_par_ba += two_op.parameters(Spin::beta, Spin::alpha).Eigen();
                total_two_op_par_bb += two_op.parameters(Spin::beta, Spin::beta).Eigen();
            }
            this->total_two_op = ScalarUSQTwoElectronOperator<Scalar>(total_two_op_par_aa, total_two_op_par_ab, total_two_op_par_ba, total_two_op_par_bb);
            break;
        }
        }
    }


    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const std::vector<ScalarSQTwoElectronOperator<Scalar>>& twoElectronContributionsMixed() const { return this->two_ops.allParameters(Spin::alpha, Spin::beta) + this->two_ops.allParameters(Spin::beta, Spin::alpha); }

    /**
     *  @return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar>& twoElectronMixed() const { return this->total_two_op.parameters(Spin::alpha, Spin::beta) + this->total_two_op.parameters(Spin::beta, Spin::alpha); }
};  // namespace GQCP


}  // namespace GQCP
