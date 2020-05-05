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


#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Basis/SpinorBasis/Spin.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Utilities/typedefs.hpp"


namespace GQCP {


/**
 *  A class that represents a spinor basis without any restrictions (G for generalized) on the expansion of the alpha and beta components in terms of the underlying (possibly different) scalar bases
 * 
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell the underlying scalar bases contain
 */
template <typename _ExpansionScalar, typename _Shell>
class GSpinorBasis: public SimpleSpinorBasis<_ExpansionScalar, GSpinorBasis<_ExpansionScalar, _Shell>> {
public:
    using ExpansionScalar = _ExpansionScalar;
    using Shell = _Shell;

    using Base = SimpleSpinorBasis<_ExpansionScalar, GSpinorBasis<_ExpansionScalar, _Shell>>;


private:
    std::array<ScalarBasis<Shell>, 2> scalar_bases;  // the scalar bases for the alpha and beta components


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param alpha_scalar_basis           the scalar basis in which the alpha components are expanded
     *  @param beta_scalar_basis            the scalar basis in which the beta components are expanded
     *  @param C                            the coefficient matrix, i.e. the matrix of the expansion coefficients of the spinors in terms of the underlying scalar basis
     */
    GSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        Base(C),
        scalar_bases {alpha_scalar_basis, beta_scalar_basis} {
        // Check if the dimensions of the given objects are compatible
        const auto K_alpha = alpha_scalar_basis.numberOfBasisFunctions();
        const auto K_beta = beta_scalar_basis.numberOfBasisFunctions();

        if (C.dimension() != K_alpha + K_beta) {
            throw std::invalid_argument("GSpinorBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const TransformationMatrix<ExpansionScalar>&): The given dimensions of the scalar bases and coefficient matrix are incompatible.");
        }
    }


    /**
     *  Construct a generalized spinor basis in which both underlying scalar bases are equal
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     *  @param C                    the coefficient matrix, i.e. the matrix of the expansion coefficients of the spinors in terms of the underlying scalar basis
     */
    GSpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        GSpinorBasis(scalar_basis, scalar_basis, C) {}


    /**
     *  Construct a generalized spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity
     */
    GSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        GSpinorBasis(alpha_scalar_basis, beta_scalar_basis,
                     TransformationMatrix<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions() + beta_scalar_basis.numberOfBasisFunctions(),
                                                                     alpha_scalar_basis.numberOfBasisFunctions() + beta_scalar_basis.numberOfBasisFunctions())) {}


    /**
     *  Construct a generalized spinor basis in which both underlying scalar bases are equal, and a coefficient matrix being the identity
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    GSpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        GSpinorBasis(scalar_basis, scalar_basis) {}


    /**
     *  Construct a generalized spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    GSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct a generalized spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *
     *  @param molecule                 the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name)) {}


    /**
     *  Construct a generalized spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework            the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name_alpha          the name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component
     *  @param basisset_name_beta           the name of the basisset, e.g. "STO-3G", used for the expansion of the beta component
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    GSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha),
                     ScalarBasis<Shell>(nuclear_framework, basisset_name_beta)) {}


    /**
     *  Construct a generalized spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the molecule
     *
     *  @param molecule                     the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name_alpha          the name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component
     *  @param basisset_name_beta           the name of the basisset, e.g. "STO-3G", used for the expansion of the beta component
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha),
                     ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta)) {}


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  Convert a restricted spinor basis into a generalized framework.
     * 
     *  @param r_spinor_basis           the restricted spinor basis
     * 
     *  @return the restricted spinor basis as a generalized one
     */
    static GSpinorBasis<ExpansionScalar, Shell> FromRestricted(const RSpinorBasis<ExpansionScalar, Shell>& r_spinor_basis) {

        // Build up the 'general' coefficient matrix.
        const auto K = r_spinor_basis.numberOfSpatialOrbitals();
        const auto M = r_spinor_basis.numberOfSpinors();
        TransformationMatrix<ExpansionScalar> C_general = TransformationMatrix<ExpansionScalar>::Zero(M, M);

        C_general.topLeftCorner(K, K) = r_spinor_basis.coefficientMatrix();
        C_general.bottomRightCorner(K, K) = r_spinor_basis.coefficientMatrix();

        return GSpinorBasis<ExpansionScalar, Shell>(r_spinor_basis.scalarBasis(), C_general);  // the alpha- and beta- scalar bases are equal
    }


    /*
     *  PUBLIC METHODS
     */

    using Base::coefficientMatrix;

    /**
     *  @param sigma        alpha or beta
     * 
     *  @return the coefficient matrix for the requested spin component, i.e. the matrix of the expansion coefficients of the requested components of the spinors in terms of its underlying scalar basis
     */
    MatrixX<ExpansionScalar> coefficientMatrix(Spin sigma) const {

        const size_t K = this->numberOfCoefficients(sigma);

        switch (sigma) {
        case Spin::alpha: {
            return this->coefficientMatrix().topRows(K);
        }

        case Spin::beta: {
            return this->coefficientMatrix().bottomRows(K);
        }
        }
    }

    /**
     *  @param sigma        alpha or beta
     * 
     *  @return the scalar basis that is used for the expansion of the given spin component
     */
    const ScalarBasis<Shell>& scalarBasis(const Spin& sigma) const { return this->scalar_bases[sigma]; }

    /**
     *  @param sigma        alpha or beta
     * 
    *  @return the number of coefficients that are used for the expansion of the requested spin-component of a spinor
     */
    size_t numberOfCoefficients(const Spin& sigma) const { return this->scalarBasis(sigma).numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this generalized spinor basis
     */
    size_t numberOfSpinors() const {

        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);

        return K_alpha + K_beta;
    }

    /**
     *  @param fq_one_op        a spin-independent first-quantized operator (i.e. whose two-component matrix operator form contains the same scalar operator in the upper-left and lower-right corners)
     * 
     *  @return the second-quantized operator corresponding to the given spin-independent first-quantized operator
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> SQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, FQOneElectronOperator::Components> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, FQOneElectronOperator::Components>;


        // The strategy for calculating the matrix representation of the one-electron operator in this spinor basis is to:
        //  1. Express the operator in the underlying scalar bases; and
        //  2. Afterwards transform them using the current coefficient matrix.
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();
        QCMatrix<ResultScalar> f = QCMatrix<ResultScalar>::Zero(M, M);  // the total result

        // 1. Express the operator in the underlying scalar bases: spin-independent operators only have alpha-alpha and beta-beta blocks.
        const auto F_alpha = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis(Spin::alpha));
        const auto F_beta = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis(Spin::beta));

        f.topLeftCorner(K_alpha, K_alpha) = F_alpha;
        f.bottomRightCorner(K_beta, K_beta) = F_beta;

        // 2. Transform using the current coefficient matrix.
        ResultOperator op {{f}};  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }


    /**
     *  @param fq_one_op        the (first-quantized) electronic spin operator
     * 
     *  @return the electronic spin operator expressed in this spinor basis
     */
    auto quantize(const ElectronicSpinOperator& fq_one_op) const -> SQOneElectronOperator<product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>, ElectronicSpinOperator::Components> {

        using ResultScalar = product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, ElectronicSpinOperator::Components>;

        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();

        // The strategy to quantize the spin operator is as follows.
        //  1. First, calculate the necessary overlap integrals over the scalar bases.
        //  2. Then, construct the scalar basis representations of the components of the spin operator by placing the overlaps into the correct blocks.
        //  3. Transform the components (in scalar basis) with the current coefficient matrix to yield the components in spinor basis.

        QCMatrix<ResultScalar> S_x = QCMatrix<ResultScalar>::Zero(M, M);
        QCMatrix<ResultScalar> S_y = QCMatrix<ResultScalar>::Zero(M, M);
        QCMatrix<ResultScalar> S_z = QCMatrix<ResultScalar>::Zero(M, M);


        // 1. Calculate the necessary overlap integrals over the scalar bases.
        const auto S_aa = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalarBasis(Spin::alpha));
        const auto S_ab = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalarBasis(Spin::alpha), this->scalarBasis(Spin::beta));
        const auto S_ba = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalarBasis(Spin::beta), this->scalarBasis(Spin::alpha));
        const auto S_bb = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalarBasis(Spin::beta));


        // 2. Place the overlaps into the correct blocks.
        S_x.block(0, K_alpha, K_alpha, K_beta) = 0.5 * S_ab;
        S_x.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * S_ba;

        const cd ii(0.0, 1.0);  // 'cd' is a typedef for 'std::complex<double>', so 'ii' is the imaginary unit
        S_y.block(0, K_alpha, K_alpha, K_beta) = -0.5 * ii * S_ab;
        S_y.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * ii * S_ba;

        S_z.topLeftCorner(K_alpha, K_alpha) = 0.5 * S_aa;
        S_z.bottomRightCorner(K_beta, K_beta) = -0.5 * S_bb;


        // 3. Transform using the coefficient matrix
        ResultOperator spin_op {std::array<QCMatrix<ResultScalar>, 3> {S_x, S_y, S_z}};  // 'op' for operator
        spin_op.transform(this->coefficientMatrix());
        return spin_op;
    }


    /**
     *  @param fq_two_op        the first-quantized Coulomb operator
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_two_op) const -> SQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Components> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Components>;

        // The strategy for calculating the matrix representation of the two-electron operator in this spinor basis is to:
        //  1. Calculate the Coulomb integrals in the underlying scalar bases;
        //  2. Place the calculated integrals as 'blocks' in the larger representation, so that we can;
        //  3. Transform the operator using the current coefficient matrix.

        // 1. Calculate the Coulomb integrals in the underlying scalar bases.
        const auto g_aaaa = IntegralCalculator::calculateLibintIntegrals(Operator::Coulomb(), this->scalarBasis(Spin::alpha));
        const auto g_aabb = IntegralCalculator::calculateLibintIntegrals(Operator::Coulomb(), this->scalarBasis(Spin::alpha), this->scalarBasis(Spin::beta));
        const auto g_bbaa = IntegralCalculator::calculateLibintIntegrals(Operator::Coulomb(), this->scalarBasis(Spin::beta), this->scalarBasis(Spin::alpha));
        const auto g_bbbb = IntegralCalculator::calculateLibintIntegrals(Operator::Coulomb(), this->scalarBasis(Spin::beta));


        // 2. Place the calculated integrals as 'blocks' in the larger representation
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);

        const auto M = this->numberOfSpinors();
        QCRankFourTensor<ResultScalar> g_par(M);  // 'par' for 'parameters'
        g_par.setZero();

        // Primed indices are indices in the larger representation, normal ones are those in the smaller tensors.
        for (size_t mu_ = 0; mu_ < M; mu_++) {  // mu 'prime'
            const size_t mu = mu_ % K_alpha;

            for (size_t nu_ = 0; nu_ < M; nu_++) {  // nu 'prime'
                const size_t nu = nu_ % K_alpha;

                for (size_t rho_ = 0; rho_ < M; rho_++) {  // rho 'prime'
                    const size_t rho = rho_ % K_alpha;

                    for (size_t lambda_ = 0; lambda_ < M; lambda_++) {  // lambda 'prime'
                        const size_t lambda = lambda_ % K_alpha;

                        if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_aaaa(mu, nu, rho, lambda);
                        } else if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_aabb(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_bbaa(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_bbbb(mu, nu, rho, lambda);
                        }
                    }
                }
            }
        }


        // 3. Transform the operator using the current coefficient matrix.
        ResultOperator g_op {g_par};  // 'op' for 'operator'
        g_op.transform(this->coefficientMatrix());
        return g_op;
    }
};


}  // namespace GQCP
