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


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"


namespace GQCP {


/**
 *  A class that represents an unrestricted spinor basis (U for unrestricted), the alpha and beta components have an individual (possibly different) expansion in their (possibly different) underlying scalar bases
 * 
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell the underlying scalar bases contain
 */
template <typename _ExpansionScalar, typename _Shell>
class USpinorBasis {
public:
    using ExpansionScalar = _ExpansionScalar;
    using Shell = _Shell;


private:
    std::array<RSpinorBasis<ExpansionScalar, Shell>, 2> spinor_bases;  // array that holds the individual spinor bases for the alpha and beta components (in that order)


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param alpha_scalar_basis           the scalar basis in which the alpha components are expanded
     *  @param beta_scalar_basis            the scalar basis in which the beta components are expanded
     *  @param C_alpha                      the alpha coefficient matrix, i.e. the matrix of the expansion coefficients of the alpha spinors in terms of the underlying scalar basis
     *  @param C_beta                       the beta coefficient matrix, i.e. the matrix of the expansion coefficients of the beta spinors in terms of the underlying scalar basis
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const TransformationMatrix<ExpansionScalar>& C_alpha, const TransformationMatrix<ExpansionScalar>& C_beta) :
        spinor_bases ({RSpinorBasis<ExpansionScalar, Shell>(alpha_scalar_basis, C_alpha), RSpinorBasis<ExpansionScalar, Shell>(beta_scalar_basis, C_beta)})
    {
        // Check if the dimensions of the given objects are compatible
        const auto K_alpha = alpha_scalar_basis.numberOfBasisFunctions();
        const auto K_beta = beta_scalar_basis.numberOfBasisFunctions();

        if (C_alpha.dimension() != K_alpha) {
            throw std::invalid_argument("USpinorBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const TransformationMatrix<ExpansionScalar>&, const TransformationMatrix<ExpansionScalar>): The given dimensions of the scalar basis and coefficient matrix for alpha are incompatible.");
        }

        if (C_beta.dimension() != K_beta) {
            throw std::invalid_argument("USpinorBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const TransformationMatrix<ExpansionScalar>&, const TransformationMatrix<ExpansionScalar>): The given dimensions of the scalar basis and coefficient matrix for beta are incompatible.");
        }
    }


    /**
     *  Construct an unrestricted spinor basis in which both underlying scalar bases and their expansions are equal
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     *  @param C                    the coefficient matrix, i.e. the matrix of the expansion coefficients of the spinors in terms of the underlying scalar bases
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        USpinorBasis(scalar_basis, scalar_basis, C, C)
    {}


    /**
     *  Construct an unrestricted spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) : 
        USpinorBasis(alpha_scalar_basis, beta_scalar_basis, TransformationMatrix<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions(), alpha_scalar_basis.numberOfBasisFunctions()), TransformationMatrix<ExpansionScalar>::Identity(beta_scalar_basis.numberOfBasisFunctions(), beta_scalar_basis.numberOfBasisFunctions()))
    {}


    /**
     *  Construct an unrestricted spinor basis in which both underlying scalar bases are equal, and the coefficient matrices are the identity
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        USpinorBasis(scalar_basis, scalar_basis)
    {}


    /**
     *  Construct an unrestricted spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *      
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting unrestricted spinor basis is (most likely) non-orthogonal
     * 
     */
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name))
    {}


    /**
     *  Construct an unrestricted spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *      
     *  @param molecule                 the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name))
    {}


    /**
     *  Construct an unrestricted spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework            the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name_alpha          the name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component
     *  @param basisset_name_beta           the name of the basisset, e.g. "STO-3G", used for the expansion of the beta component
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting unrestricted spinor basis is (most likely) non-orthogonal
     */
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha), ScalarBasis<Shell>(nuclear_framework, basisset_name_beta))
    {}


    /**
     *  Construct an unrestricted spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the molecule
     *
     *  @param molecule                     the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name_alpha          the name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component
     *  @param basisset_name_beta           the name of the basisset, e.g. "STO-3G", used for the expansion of the beta component
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting unrestricted spinor basis is (most likely) non-orthogonal
     */
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha), ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta))
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param component        the spin component
     * 
     *  @return the coefficient matrix for the requested component, i.e. the matrix of the expansion coefficients of the requested component of the spinors in terms of its underlying scalar basis
     */
    const MatrixX<ExpansionScalar>& coefficientMatrix(SpinComponent component) const { 
        return spinor_bases[component].coefficientMatrix();
    }

    /**
     *  @param component        the spin component
     * 
     *  @return the scalar basis in which the requested component is expanded
     */
    const ScalarBasis<Shell>& scalarBasis(const SpinComponent& component) const { return this->spinor_bases[component].scalarBasis(); }

    /**
     *  @param component        the spin component
     * 
     *  @return the underlying spinor basis for a single component
     */
    const RSpinorBasis<ExpansionScalar, Shell>& spinorBasis(const SpinComponent& component) const { return this->spinor_bases[component]; }

    /**
     *  @param component        the spin component
     * 
     *  @return the scalar basis in which the requested component is expanded
     */
    size_t numberOfCoefficients(const SpinComponent& component) const { return this->scalarBasis(component).numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this unrestricted spinor basis
     */
    size_t numberOfSpinors() const { 

        const auto K_alpha = this->numberOfCoefficients(SpinComponent::ALPHA);
        const auto K_beta = this->numberOfCoefficients(SpinComponent::BETA);

        return K_alpha + K_beta;
    }

    /**
     *  @param component                the spin component
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this spinor basis for the requested component is orthonormal within the given precision
     */
    bool isOrthonormal(const SpinComponent& component, const double precision = 1.0e-08) const {
        return this->spinor_bases[component].isOrthonormal();
    }

    /**
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this spinor basis is orthonormal within the given precision
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {
        return this->isOrthonormal(SpinComponent::ALPHA, precision) && this->isOrthonormal(SpinComponent::BETA, precision);
    }

    /**
     *  @param component                the spin component
     * 
     *  @return the transformation matrix to the Löwdin basis for the requested component: T = S_current^{-1/2}
     */
    TransformationMatrix<double> lowdinOrthonormalizationMatrix(const SpinComponent& component) const {
        return this->scalarBasis(component).lowdinOrthonormalizationMatrix();
    }

    /**
     *  @param component                the spin component
     * 
     *  Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix
     */
    void lowdinOrthonormalize(const SpinComponent& component) {
        this->spinor_bases[component].lowdinOrthonormalize();
    }

    /**
     *  Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix
     */
    void lowdinOrthonormalize() {
        this->spinor_bases[SpinComponent::ALPHA].lowdinOrthonormalize();
        this->spinor_bases[SpinComponent::BETA].lowdinOrthonormalize();
    }

    /**
     *  @param component                the spin component
     * 
     *  @return the overlap (one-electron) operator of the requested component of this spinor basis
     */
    ScalarSQOneElectronOperator<ExpansionScalar> overlap(const SpinComponent& component) const {
        return this->spinor_bases[component].quantize(Operator::Overlap());
    }


    /**
     *  Rotate the spinor basis of the requested component to another one using the given unitary transformation matrix
     * 
     *  @param U                        the unitary transformation matrix that transforms both the alpha- and beta components
     *  @param component                the spin component
     */
    void rotate(const TransformationMatrix<ExpansionScalar>& U, const SpinComponent& component) {
        this->spinor_bases[component].rotate();
    }


    /**
     *  Rotate the spinor basis to another one using the given unitary transformation matrix
     * 
     *  @param U            the unitary transformation matrix that transforms both the alpha- and beta components
     */
    void rotate(const TransformationMatrix<ExpansionScalar>& U) {
        this->rotate(SpinComponent::ALPHA);
        this->rotate(SpinComponent::BETA);
    }


    /**
     *  Rotate the spinor basis of the requested component to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     *  @param component                the spin component
     * 
     *  @note this function is only available for real spinor bases because Jacobi rotation parameters generate real rotations
     */
    template<typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters, const SpinComponent& component) {
        this->spinor_bases[component].rotate(jacobi_rotation_parameters);
    }


    /**
     *  Rotate the spinor basis for both the alpha- and beta components to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     * 
     *  @note this function is only available for real spinor bases because Jacobi rotation parameters generate real rotations
     */
    template<typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {
        this->rotate(jacobi_rotation_parameters, SpinComponent::ALPHA);
        this->rotate(jacobi_rotation_parameters, SpinComponent::BETA);
    }


    /**
     *  Transform the spinor basis for one component to another one using the given transformation matrix
     *
     *  @param T                        the transformation matrix that transforms the requested component
     *  @param component                the spin component
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T, const SpinComponent& component) {
         this->spinor_bases[component].transform(T);
    }


    /**
     *  Transform the spinor basis to another one using the given transformation matrix
     *
     *  @param T            the transformation matrix that transforms both the alpha- and beta components
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T) {
         this->transform(T, SpinComponent::ALPHA);
         this->transform(T, SpinComponent::BETA);
    }


    /**
     *  @param ao_list          indices of the AOs used for the Mulliken populations
     *  @param component        the spin component
     *
     *  @return the Mulliken operator for a set of AOs and the requested component
     *
     *  @note this method is only available for real matrix representations
     */
    template<typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const Vectoru& ao_list, const SpinComponent& component) const {
        return this->spinor_bases[component].template calculateMullikenOperator<ExpansionScalar>(ao_list);
    }


    /**
     *  @param fq_op            the first-quantized one-electron operator
     *  @param component        the spin component
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator in the spinor basis of the requested component
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_op, const SpinComponent& component) const -> SQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, FQOneElectronOperator::Components> {
        return this->spinor_bases[component].quantize(fq_op);
    }
        

    /**
     *  @param fq_op            the first-quantized Coulomb operator
     *  @param component        the spin component
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator in the spinor basis of the requested component
     */
    auto quantize(const CoulombRepulsionOperator& fq_op, const SpinComponent& component) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Components> {
        return this->spinor_bases[component].quantize(fq_op);
    }
};


}  // namespace GQCP
