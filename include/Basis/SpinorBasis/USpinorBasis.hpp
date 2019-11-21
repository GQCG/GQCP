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


#include "Basis/SpinorBasis/CompoundSpinorBasis.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"



namespace GQCP {


/**
 *  A class that represents an unrestricted spinor basis (U for unrestricted), the alpha and beta components have an individual (possibly different) expansion in their (possibly different) underlying scalar bases
 * 
 *  @tparam _ExpansionScalar        the scalar type of the expansion coefficients
 *  @tparam _Shell                  the type of shell the underlying scalar bases contain
 */
template <typename _ExpansionScalar, typename _Shell>
class USpinorBasis : public CompoundSpinorBasis<_Shell> {
private:
    std::array<TransformationMatrix<ExpansionScalar>, 2> C_array;  // array that holds the expansion coefficients for the alpha and beta components


public:
    using ExpansionScalar = _ExpansionScalar;
    using Shell = _Shell;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param alpha_scalar_basis           the scalar basis in which the alpha components are expanded
     *  @param beta_scalar_basis            the scalar basis in which the beta components are expanded
     *  @param C_alpha                      the coefficient matrix, i.e. the matrix of the expansion coefficients of the alpha spinors in terms of the underlying scalar basis
     *  @param C_beta                       the coefficient matrix, i.e. the matrix of the expansion coefficients of the beta spinors in terms of the underlying scalar basis
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const TransformationMatrix<ExpansionScalar>& C_alpha, const TransformationMatrix<ExpansionScalar>& C_beta) :
        CompoundSpinorBasis(alpha_scalar_basis, beta_scalar_basis),
        C_array ({C_alpha, C_beta})
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
     *  Construct a unrestricted spinor basis in which both underlying scalar bases and their expansions are equal
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     *  @param C                    the coefficient matrix, i.e. the matrix of the expansion coefficients of the spinors in terms of the underlying scalar basis
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis, const TransformationMatrix<ExpansionScalar>& C) :
        USpinorBasis(scalar_basis, scalar_basis, C, C)
    {}


    /**
     *  Construct a unrestricted spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) : 
        USpinorBasis(alpha_scalar_basis, beta_scalar_basis, TransformationMatrix<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions(), alpha_scalar_basis.numberOfBasisFunctions()), TransformationMatrix<ExpansionScalar>::Identity(beta_scalar_basis.numberOfBasisFunctions(), beta_scalar_basis.numberOfBasisFunctions())
    {}


    /**
     *  Construct a unrestricted spinor basis in which both underlying scalar bases are equal, and the coefficient matrices are the identity
     * 
     *  @param scalar_basis         the scalar basis in which both the alpha and beta components are expanded
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        USpinorBasis(scalar_basis, scalar_basis)
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
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name))
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
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name))
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
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha), ScalarBasis<Shell>(nuclear_framework, basisset_name_beta))
    {}


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
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha), ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta))
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param component        the spin component
     * 
     *  @return the coefficient matrix for the requested component, i.e. the matrix of the expansion coefficients of the requested components of the spinors in terms of its underlying scalar basis
     */
    const MatrixX<ExpansionScalar>& coefficientMatrix(SpinComponent component) const { 
        return C_array[component];
    }
};


}  // namespace GQCP
