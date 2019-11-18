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
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"


namespace GQCP {


/**
 *  A class that represents a spinor basis without any restrictions (G for generalized) on the expansion of the alpha and beta components in terms of the underlying (possibly different) scalar bases
 * 
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell the underlying scalar bases contain
 */
template <typename _ExpansionScalar, typename _Shell>
class GSpinorBasis : public SimpleSpinorBasis<_ExpansionScalar, GSpinorBasis<_ExpansionScalar, _Shell>> {
public:
    using ExpansionScalar = _ExpansionScalar;
    using Shell = _Shell;


private:
    ScalarBasis<Shell> alpha_scalar_basis;  // the scalar basis in which the alpha components are expanded
    ScalarBasis<Shell> beta_scalar_basis;  // the scalar basis in which the beta components are expanded


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
        SimpleSpinorBasis<ExpansionScalar>(C),
        alpha_scalar_basis (alpha_scalar_basis),
        beta_scalar_basis (beta_scalar_basis)
    {
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
        GSpinorBasis(scalar_basis, scalar_basis, C)
    {}


    /**
     *  Construct a generalized spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity
     */
    GSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) : 
        GSpinorBasis(alpha_scalar_basis, beta_scalar_basis, TransformationMatrix<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions() + beta_scalar_basis.numberOfBasisFunctions(), alpha_scalar_basis.numberOfBasisFunctions() + beta_scalar_basis.numberOfBasisFunctions()))
    {}


    /**
     *  Construct a generalized spinor basis in which both underlying scalar bases are equal, and a coefficient matrix being the identity
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    GSpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        GSpinorBasis(scalar_basis, scalar_basis)
    {}


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
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name))
    {}


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
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name))
    {}


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
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha), ScalarBasis<Shell>(nuclear_framework, basisset_name_beta))
    {}


    /**
     *  Construct a generalized spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the molecule
     *
     *  @param molecule                 the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name_alpha          the name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component
     *  @param basisset_name_beta           the name of the basisset, e.g. "STO-3G", used for the expansion of the beta component
     *
     *  @note the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     *  @note the resulting generalized spinor basis is (most likely) non-orthogonal
     */
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha), ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta))
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the alpha coefficient matrix, i.e. the matrix of the expansion coefficients of the alpha-components of the spinors in terms of the underlying alpha-scalar basis
     */
    MatrixX<ExpansionScalar> alphaCoefficientMatrix() const { return this->C.topRows(this->numberOfAlphaCoefficients()); }

    /**
     *  @return the scalar basis in which the alpha components are expanded
     */
    const ScalarBasis<Shell>& alphaScalarBasis() const { return this->alpha_scalar_basis; }

    /**
     *  @return the scalar basis in which the beta components are expanded
     */
    const ScalarBasis<Shell>& betaScalarBasis() const { return this->beta_scalar_basis; }

    /**
     *  @return the alpha coefficient matrix, i.e. the matrix of the expansion coefficients of the alpha-components of the spinors in terms of the underlying alpha-scalar basis
     */
    MatrixX<ExpansionScalar> betaCoefficientMatrix() const { return this->C.bottomRows(this->numberOfBetaCoefficients()); }

    /**
     *  @return the number of coefficients that are used for the expansion of the alpha-component of a spinor
     */
    size_t numberOfAlphaCoefficients() const { return this->alphaScalarBasis().numberOfBasisFunctions(); }

    /**
     *  @return the number of coefficients that are used for the expansion of the beta-component of a spinor
     */
    size_t numberOfBetaCoefficients() const { return this->betaScalarBasis().numberOfBasisFunctions(); }

    /**
     *  @return the number of spinors that 'are' in this generalized spinor basis
     */
    size_t numberOfSpinors() const { 
        
        const auto K_alpha = this->numberOfAlphaCoefficients();
        const auto K_beta = this->numberOfBetaCoefficients();

        return K_alpha + K_beta;
    }


    /**
     *  @param fq_op        the first-quantized operator
     * 
     *  @return the second-quantized operator corresponding to the given first-quantized operator
     */
    auto quantize(const OverlapOperator& fq_op) const -> SQOneElectronOperator<product_t<OverlapOperator::Scalar, ExpansionScalar>, OverlapOperator::Components> {

        using ResultScalar = product_t<OverlapOperator::Scalar, ExpansionScalar>;
        using ResultOperator = SQOneElectronOperator<ResultScalar, OverlapOperator::Components>;


        // The strategy for calculating the matrix representation of the operator in this spinor basis is to express the operator in the underlying scalar bases and afterwards transform them using the current coefficient matrix
        const auto K_alpha = this->numberOfAlphaCoefficients();
        const auto K_beta = this->numberOfBetaCoefficients();
        const auto M = this->numberOfSpinors();
        QCMatrix<ResultScalar> f = QCMatrix<ResultScalar>::Zero(M, M);

        // Express the operator in the underlying bases: the overlap operator only has alpha-alpha and beta-beta blocks
        const auto S_alpha = this->alphaScalarBasis().calculateLibintIntegrals(fq_op);
        const auto S_beta = this->betaScalarBasis().calculateLibintIntegrals(fq_op);

        f.topLeftCorner(K_alpha, K_alpha) = S_alpha;
        f.bottomRightCorner(K_beta, K_beta) = S_beta;

        // Transform using the current coefficient matrix
        ResultOperator op ({f});  // op for 'operator'
        op.transform(this->coefficientMatrix());
        return op;
    }
};


}  // namespace GQCP
