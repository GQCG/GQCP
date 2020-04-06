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


#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/TwoRDM.hpp"
#include "Utilities/miscellaneous.hpp"
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
     *  @param one_ops       the unrestricted core (i.e. one electron) contributions to the Hamiltonian
     *  @param two_ops       the unrestricted two electron contributions to the Hamiltonian
     */
    USQHamiltonian(const std::vector<ScalarUSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarUSQTwoElectronOperator<Scalar>>& two_ops) :
        one_ops (one_ops),
        two_ops (two_ops)
    {
        // Check if the dimensions are compatible
        const auto dim = one_ops[0].dimension(GQCP::SpinComponent::ALPHA);

        if (one_ops[0].dimension(GQCP::SpinComponent::BETA) != dim) {
            throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const std::vector<ScalarUSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarUSQTwoElectronOperator<Scalar>>& two_ops: The dimensions of the alpha and beta Hamiltonian are incompatible");
        }
        
        for (const auto& two_op : this->two_ops) {
            if (two_op.dimension(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA) != dim) {
                throw std::invalid_argument("USQHamiltonian::USQHamiltonian(const std::vector<ScalarUSQOneElectronOperator<Scalar>>& one_ops, const std::vector<ScalarUSQTwoElectronOperator<Scalar>>& two_ops): The dimensions of the two electron operators are incompatible with the Hamiltonian");
            }
        }


        // Calulate the total one-electron operator
        QCMatrix<Scalar> total_one_op_par_a (dim);
        QCMatrix<Scalar> total_one_op_par_b (dim);
        total_one_op_par_a.setZero();
        total_one_op_par_b.setZero();
        for (const auto& one_op : this->one_ops) {
            total_one_op_par_a += one_op.parameters(GQCP::SpinComponent::ALPHA);
            total_one_op_par_b += one_op.parameters(GQCP::SpinComponent::BETA);
        }
        this->total_one_op = ScalarUSQOneElectronOperator<Scalar> (total_one_op_par_a, total_one_op_par_b);
        // Calculate the total two-electron operator
        QCRankFourTensor<Scalar> total_two_op_par_aa (dim);
        QCRankFourTensor<Scalar> total_two_op_par_ab (dim);
        QCRankFourTensor<Scalar> total_two_op_par_ba (dim);
        QCRankFourTensor<Scalar> total_two_op_par_bb (dim);
        total_two_op_par_aa.setZero();
        total_two_op_par_ab.setZero();
        total_two_op_par_ba.setZero();
        total_two_op_par_bb.setZero();

        for (const auto& two_op : this->two_ops) {
            total_two_op_par_aa += two_op.parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::ALPHA).Eigen();
            total_two_op_par_ab += two_op.parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA).Eigen();
            total_two_op_par_ba += two_op.parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::ALPHA).Eigen();
            total_two_op_par_bb += two_op.parameters(GQCP::SpinComponent::BETA, GQCP::SpinComponent::BETA).Eigen();
        }
        this->total_two_op = ScalarUSQTwoElectronOperator<Scalar> (total_two_op_par_aa, total_two_op_par_ab, total_two_op_par_ba, total_two_op_par_bb);
    }

     /**
     *  @param h            the (total) one-electron (i.e. core) integrals
     *  @param g            the (total) two-electron integrals
     */
    USQHamiltonian(const ScalarUSQOneElectronOperator<Scalar>& h, const ScalarUSQTwoElectronOperator<Scalar>& g) :
        USQHamiltonian(std::vector<ScalarUSQOneElectronOperator<Scalar>>{h}, std::vector<ScalarUSQTwoElectronOperator<Scalar>>{g})
    {}


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
    static enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> Molecular(const USpinorBasis<Z, GTOShell>& spinor_basis, const Molecule& molecule) {

        // Calculate the integrals for the molecular Hamiltonian
        const auto T_a = spinor_basis.quantize(Operator::Kinetic(), GQCP::SpinComponent::ALPHA);
        const auto T_b = spinor_basis.quantize(Operator::Kinetic(), GQCP::SpinComponent::BETA);

        const auto V_a = spinor_basis.quantize(Operator::NuclearAttraction(molecule), GQCP::SpinComponent::ALPHA);
        const auto V_b = spinor_basis.quantize(Operator::NuclearAttraction(molecule), GQCP::SpinComponent::BETA);

        ScalarSQOneElectronOperator<double> H_a = T_a + V_a;
        ScalarSQOneElectronOperator<double> H_b = T_b + V_b;

        GQCP::ScalarUSQOneElectronOperator<double> H (H_a.parameters(), H_b.parameters());

        const auto g_aa = spinor_basis.quantize(Operator::Coulomb(), GQCP::SpinComponent::ALPHA);
        const auto g_bb = spinor_basis.quantize(Operator::Coulomb(), GQCP::SpinComponent::BETA); 
        GQCP::QCRankFourTensor<double> g_mix (g_aa.dimension());

        for (const auto& g_values : g_aa.allParameters()) {
            g_mix += g_values.Eigen();
        }

        GQCP::ScalarUSQTwoElectronOperator<double> g (g_aa.parameters(), g_mix, g_mix, g_bb.parameters());

        return USQHamiltonian(H, g);
    }

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the dimension of the Hamiltonian, i.e. the number of spinors in which it is expressed
     */
    size_t dimension() const { return this->one_ops[0].dimension(GQCP::SpinComponent::ALPHA) + this->one_ops[0].dimension(GQCP::SpinComponent::BETA); }

    /**
     *  @return if the alpha and beta components of the unrestricted Hamiltonian are of the same dimension
     */
    bool areSpinHamiltoniansOfSameDimension() const { return this->spinHamiltonian(SpinComponent::ALPHA).dimension() == this->spinHamiltonian(SpinComponent::BETA).dimension(); }
   
    /**
     *  @param component                    the spin component
     * 
     *  @return the pure contributions of the requested component of the unrestricted Hamiltonian 
     */
    const  SQHamiltonian<Scalar> spinHamiltonian(const SpinComponent& component) const { 

        QCMatrix<Scalar> one_e = this->total_one_op.parameters(component);
        QCRankFourTensor<Scalar> two_e = this->total_two_op.parameters(component, component);

        return SQHamiltonian<Scalar>(one_e, two_e); 
    }

    /**
     *  @return the total mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const ScalarSQTwoElectronOperator<Scalar> twoElectronMixed() const { 
        QCRankFourTensor<Scalar> two_op_mix = this->total_two_op.parameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA); 
        return ScalarSQTwoElectronOperator<Scalar>(two_op_mix);
    }

    /**
     *  @return the contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian
     */
    const std::vector<ScalarSQTwoElectronOperator<Scalar>>& twoElectronContributionsMixed() const { 
        return this->two_ops.allParameters(GQCP::SpinComponent::ALPHA, GQCP::SpinComponent::BETA); 
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
     *  @param component                the spin component
     */
    void transform(const TransformationMatrix<Scalar>& T, const SpinComponent& component) {

        // Transform the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.parameters(component).basisTransformInPlace(T);
        }

        // transform the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.parameters(component, component).basisTransformInPlace(T);
        }

        // Transform the totals
        this->total_one_op.parameters(component).basisTransformInPlace(T);
        this->total_two_op.parameters(component, component).basisTransformInPlace(T);
    }


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
     *  @param component            the spin component
     */
    void rotate(const TransformationMatrix<Scalar>& U, const SpinComponent& component) {

        // Rotate the one-electron contributions
        for (auto& one_op : this->one_ops) {
            one_op.allParameters(component).basisRotateInPlace(U);
        }

        // Rotate the two-electron contributions
        for (auto& two_op : this->two_ops) {
            two_op.allParameters(component, component).basisRotateInPlace(U);
        }

        // Rotate the totals
        this->total_one_op.allParameters(component).basisRotateInPlace(U);
        this->total_two_op.allParameters(component, component).basisRotateInPlace(U);
    }
    

    /**
     *  In-place rotate the matrix representations of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {
        this->rotate(jacobi_rotation_parameters);
    }


    /**
     *  In-place rotate the matrix representation of one of the spin components of the Hamiltonian using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. Note that this function is only available for real (double) matrix representations
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     *  @param component                        the spin component
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters, const SpinComponent& component) {
        this->rotate(jacobi_rotation_parameters, component);
    }


    /**
     *  Constrain a spin component of the unrestricted Hamiltonian according to the convention: - lambda * constraint
     *
     *  @param one_op           the one-electron operator used as a constraint
     *  @param lambda           the Lagrangian multiplier for the constraint
     *  @param component        the spin component
     *
     *  @return a copy of the constrained Hamiltonian
     *
     *  Note that this method is only available for real matrix representations
     */
    template<typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, USQHamiltonian<double>> constrain(const ScalarSQOneElectronOperator<double>& one_op, const double lambda, const SpinComponent& component) const {

        auto const constrained_one_ops = this->one_ops.allParameters(component) - lambda * one_op;
        return USQHamiltonian(constrained_one_ops, this->two_ops);
        
    }
};


}  // namespace GQCP
