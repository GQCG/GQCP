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


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/Spin.hpp"
#include "Utilities/type_traits.hpp"


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
        spinor_bases {RSpinorBasis<ExpansionScalar, Shell>(alpha_scalar_basis, C_alpha),
                      RSpinorBasis<ExpansionScalar, Shell>(beta_scalar_basis, C_beta)} {

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
        USpinorBasis(scalar_basis, scalar_basis, C, C) {}


    /**
     *  Construct an unrestricted spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        USpinorBasis(alpha_scalar_basis, beta_scalar_basis,
                     TransformationMatrix<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions(), alpha_scalar_basis.numberOfBasisFunctions()),
                     TransformationMatrix<ExpansionScalar>::Identity(beta_scalar_basis.numberOfBasisFunctions(), beta_scalar_basis.numberOfBasisFunctions())) {}


    /**
     *  Construct an unrestricted spinor basis in which both underlying scalar bases are equal, and the coefficient matrices are the identity
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        USpinorBasis(scalar_basis, scalar_basis) {}


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
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct an unrestricted spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *      
     *  @param molecule                 the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting unrestricted spinor basis is (most likely) non-orthogonal
     */
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name)) {}


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
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha),
                     ScalarBasis<Shell>(nuclear_framework, basisset_name_beta)) {}


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
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha),
                     ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta)) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Create an unrestricted spinor basis from a restricted spinor basis, leading to alpha- and beta- coefficient matrices that are equal.
     * 
     *  @param r_spinor_basis               the restricted spinor basis
     * 
     *  @return an unrestricted spinor basis
     */
    static USpinorBasis<ExpansionScalar, Shell> FromRestricted(const RSpinorBasis<ExpansionScalar, Shell>& r_spinor_basis) {

        const auto scalar_basis = r_spinor_basis.scalarBasis();
        const auto C = r_spinor_basis.coefficientMatrix();
        return USpinorBasis<ExpansionScalar, Shell>(scalar_basis, scalar_basis, C, C);
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param ao_list              indices of the AOs used for the atomic spin operator in the z-direction
     *  @param sigma                alpha or beta
     *
     *  @return the SQ atomic spin operator in the z-direction for a set of AOs
     *
     *  Note that this method is only available for real SQoperators
     */
    template <typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateAtomicSpinZ(const std::vector<size_t>& ao_list, const Spin& sigma) const {

        // The atomic spin operator can be calculated as as the atomic Mulliken operator divided by 2, multiplied by the correct sign factor
        int sign = 1 - 2 * sigma;  // 1 for ALPHA, -1 for BETA
        const auto spin_z_par = 0.5 * sign * this->spinor_bases[sigma].calculateMullikenOperator(ao_list).parameters();
        return ScalarSQOneElectronOperator<double>(spin_z_par);
    }


    /**
     *  @param ao_list              indices of the AOs used for the Mulliken populations
     *  @param sigma                alpha or beta
     *
     *  @return the Mulliken operator for a set of AOs and the requested component
     *
     *  @note this method is only available for real matrix representations
     */
    template <typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<double>> calculateMullikenOperator(const std::vector<size_t>& ao_list, const Spin& sigma) const {
        return this->spinor_bases[sigma].template calculateMullikenOperator<ExpansionScalar>(ao_list);
    }


    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the coefficient matrix for the requested component, i.e. the matrix of the expansion coefficients of the requested component of the spinors in terms of its underlying scalar basis
     */
    const MatrixX<ExpansionScalar>& coefficientMatrix(Spin sigma) const {
        return spinor_bases[sigma].coefficientMatrix();
    }


    /**
     *  @return the total number of spinors/spin-orbitals that this spinor basis describes
     */
    size_t numberOfSpinors() const {

        const auto K_alpha = this->numberOfSpinors(Spin::alpha);
        const auto K_beta = this->numberOfSpinors(Spin::beta);

        return K_alpha + K_beta;
    }

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the total number of sigma-spinors/spin-orbitals that this spinor basis describes
     */
    size_t numberOfSpinors(GQCP::Spin sigma) const {
        return this->scalarBasis(sigma).numberOfBasisFunctions();
    }


    /**
     *  @param sigma                alpha or beta
     *  @param precision            the precision used to test orthonormality
     * 
     *  @return if this spinor basis for the requested component is orthonormal within the given precision
     */
    bool isOrthonormal(const Spin& sigma, const double precision = 1.0e-08) const {
        return this->spinor_bases[sigma].isOrthonormal();
    }


    /**
     *  @param precision                the precision used to test orthonormality
     * 
     *  @return if this spinor basis is orthonormal within the given precision
     */
    bool isOrthonormal(const double precision = 1.0e-08) const {
        return this->isOrthonormal(Spin::alpha, precision) && this->isOrthonormal(Spin::beta, precision);
    }


    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the transformation matrix to the Löwdin basis for the requested component: T = S_current^{-1/2}
     */
    TransformationMatrix<double> lowdinOrthonormalizationMatrix(const Spin& sigma) const {
        return this->spinor_bases[sigma].lowdinOrthonormalizationMatrix();
    }


    /**
     *  @param sigma                alpha or beta
     * 
     *  Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix
     */
    void lowdinOrthonormalize(const Spin& sigma) {
        this->spinor_bases[sigma].lowdinOrthonormalize();
    }


    /**
     *  Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix
     */
    void lowdinOrthonormalize() {
        this->spinor_bases[Spin::alpha].lowdinOrthonormalize();
        this->spinor_bases[Spin::beta].lowdinOrthonormalize();
    }


    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the overlap (one-electron) operator of the requested component of this spinor basis
     */
    ScalarSQOneElectronOperator<ExpansionScalar> overlap(const Spin& sigma) const {
        return this->spinor_bases[sigma].quantize(Operator::Overlap());
    }


    /**
     *  @param fq_op                the first-quantized one-electron operator
     *  @param sigma                alpha or beta
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator in the spinor basis of the requested component
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_op, const Spin& sigma) const -> SQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, FQOneElectronOperator::Components> {
        return this->spinor_bases[sigma].quantize(fq_op);
    }


    /**
     *  @param fq_op                the first-quantized Coulomb operator
     *  @param sigma                alpha or beta
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator in the spinor basis of the requested component
     * 
     *  @note This method is not (yet) capable of calculating 'mixed' integrals such as g_aabb.
     */
    auto quantize(const CoulombRepulsionOperator& fq_op, const Spin& sigma) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Components> {
        return this->spinor_bases[sigma].quantize(fq_op);
    }


    /**
     *  Rotate the spinor basis of the requested component to another one using the given unitary transformation matrix
     * 
     *  @param U                        the unitary transformation matrix that transforms both the alpha or beta component
     *  @param sigma                    alpha or beta
     */
    void rotate(const TransformationMatrix<ExpansionScalar>& U, const Spin& sigma) {
        this->spinor_bases[sigma].rotate(U);
    }

    /**
     *  Rotate the spinor basis to another one using the given unitary transformation matrix
     * 
     *  @param U            the unitary transformation matrix that transforms both the alpha- and beta components
     * 
     *  @note this method is only valid when the beta and alpha component are of the same dimension, and will only accept matrices of the same dimension as the individual component.
     */
    void rotate(const TransformationMatrix<ExpansionScalar>& U) {
        this->rotate(Spin::alpha);
        this->rotate(Spin::beta);
    }


    /**
     *  Rotate the spinor basis of the requested component to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     *  @param sigma                            alpha or beta
     * 
     *  @note this function is only available for real spinor bases because Jacobi rotation parameters generate real rotations
     */
    template <typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters, const Spin& sigma) {
        this->spinor_bases[sigma].rotate(jacobi_rotation_parameters);
    }


    /**
     *  Rotate the spinor basis for both the alpha- and beta components to another one using the unitary transformation matrix that corresponds to the given Jacobi rotation parameters
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     * 
     *  @note this function is only available for real spinor bases because Jacobi rotation parameters generate real rotations
     *  @note this method is only valid when the beta and alpha component are of the same dimension, and will only accept parameters of the same dimension as the individual component.
     */
    template <typename Z = ExpansionScalar>
    enable_if_t<std::is_same<Z, double>::value> rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {
        this->rotate(jacobi_rotation_parameters, Spin::alpha);
        this->rotate(jacobi_rotation_parameters, Spin::beta);
    }


    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the scalar basis in which the requested component is expanded
     */
    const ScalarBasis<Shell>& scalarBasis(const Spin& sigma) const { return this->spinor_bases[sigma].scalarBasis(); }

    /**
     *  @param sigma                alpha or beta
     * 
     *  @return the underlying spinor basis for a single component
     */
    const RSpinorBasis<ExpansionScalar, Shell>& spinorBasis(const Spin& sigma) const { return this->spinor_bases[sigma]; }


    /**
     *  Transform the spinor basis for one component to another one using the given transformation matrix
     *
     *  @param T                        the transformation matrix that transforms the requested component
     *  @param sigma                    alpha or beta
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T, const Spin& sigma) {
        this->spinor_bases[sigma].transform(T);
    }


    /**
     *  Transform the spinor basis to another one using the given transformation matrices.
     *
     *  @param T_alpha              the transformation matrix that transforms the alpha- spin-orbitals
     *  @param T_beta               the transformation matrix that transforms the beta- spin-orbitals
     * 
     *  @note this method is only valid when the beta and alpha component are of the same dimension, and will only accept matrices of the same dimension as the individual component.
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T_alpha, const TransformationMatrix<ExpansionScalar>& T_beta) {
        this->transform(T_alpha, Spin::alpha);
        this->transform(T_beta, Spin::beta);
    }


    /**
     *  Transform the spinor basis to another one using the given transformation matrix
     *
     *  @param T            the transformation matrix that transforms both the alpha- and beta components
     * 
     *  @note this method is only valid when the beta and alpha component are of the same dimension, and will only accept matrices of the same dimension as the individual component.
     */
    void transform(const TransformationMatrix<ExpansionScalar>& T) {
        this->transform(T, T);
    }
};


}  // namespace GQCP
