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
#include "Basis/SpinorBasis/USpinOrbitalBasisComponent.hpp"
#include "Basis/Transformations/SpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/SpinResolvedJacobiRotatable.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  A class that represents an unrestricted spin-orbital basis. The difference with a restricted spin-orbital basis is that the alpha- and beta-spin-orbitals have an individual (i.e. possibly different) expansion in their (possibly different) underlying scalar bases.
 * 
 *  @tparam _ExpansionScalar                The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                          The type of shell that the underlying scalar basis contains.
the type of shell the underlying scalar bases contain
 */
template <typename _ExpansionScalar, typename _Shell>
class USpinorBasis:
    public SpinResolvedBase<USpinOrbitalBasisComponent<_ExpansionScalar, _Shell>, USpinorBasis<_ExpansionScalar, _Shell>>,
    public SpinResolvedBasisTransformable<USpinorBasis<_ExpansionScalar, _Shell>>,
    public SpinResolvedJacobiRotatable<USpinorBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<USpinOrbitalBasisComponent<ExpansionScalar, Shell>, USpinorBasis<ExpansionScalar, Shell>>::SpinResolvedBase;


    /**
     *  Create an `USpinorBasis` from an alpha and a beta-spin-orbital basis and a transformation that expresses the current spin-orbitals in terms of that underlying scalar basis.
     * 
     *  @param alpha_scalar_basis           The scalar basis in which the alpha spin-orbitals are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta spin-orbitals are expanded.
     *  @param C                            The transformation that expresses the current spin-orbitals in terms of the underlying scalar basis.
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const UTransformationMatrix<ExpansionScalar>& C) :
        SpinResolvedBase(USpinOrbitalBasisComponent<ExpansionScalar, Shell> {alpha_scalar_basis, C.alpha()},
                         USpinOrbitalBasisComponent<ExpansionScalar, Shell> {beta_scalar_basis, C.beta()}) {

        // Check if the dimensions of the given objects are compatible.
        if (C.alpha().numberOfOrbitals() != alpha_scalar_basis.numberOfBasisFunctions()) {
            throw std::invalid_argument("USpinorBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const UTransformationMatrix<ExpansionScalar>&): The given dimensions of the scalar basis and coefficient matrix for the alpha spin-orbitals are incompatible.");
        }

        if (C.beta().numberOfOrbitals() != beta_scalar_basis.numberOfBasisFunctions()) {
            throw std::invalid_argument("USpinorBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const UTransformationMatrix<ExpansionScalar>&): The given dimensions of the scalar basis and coefficient matrix for the beta spin-orbitals are incompatible.");
        }
    }


    /**
     *  Construct an unrestricted spin-orbital basis in which both underlying scalar bases and their expansions are equal.
     *
     *  @param scalar_basis         The scalar basis in which both the alpha and beta spin-orbitals are expanded.
     *  @param C                    The transformation that expresses the current spin-orbitals in terms of the underlying scalar basis.
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis, const UTransformationMatrixComponent<ExpansionScalar>& C) :
        USpinorBasis(scalar_basis, scalar_basis, UTransformationMatrix<ExpansionScalar>::FromEqual(C)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with two different underlying scalar bases, and a coefficient matrix being the identity. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     * 
     *  @param alpha_scalar_basis           The scalar basis in which the alpha spin-orbitals are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta spin-orbitals are expanded.
     */
    USpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        USpinorBasis(alpha_scalar_basis, beta_scalar_basis,
                     UTransformationMatrix<ExpansionScalar>(
                         UTransformationMatrixComponent<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions()),
                         UTransformationMatrixComponent<ExpansionScalar>::Identity(beta_scalar_basis.numberOfBasisFunctions()))) {}


    /**
     *  Construct an unrestricted spin-orbital basis in which both underlying scalar bases are equal, and the coefficient matrix being the identity. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param scalar_basis         The scalar basis in which both the alpha and beta spin-orbitals are expanded.
     */
    USpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        USpinorBasis(scalar_basis, scalar_basis) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases (equal for both the alpha and beta components) that are made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells of the scalar basis should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases (equal for both the alpha and beta components) that are made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param molecule                 The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the nuclear framework. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param nuclear_framework            The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name_alpha          The name of the basisset, e.g. "STO-3G", used for the expansion of the alpha spin-orbitals.
     *  @param basisset_name_beta           The name of the basisset, e.g. "STO-3G", used for the expansion of the beta spin-orbitals.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha),
                     ScalarBasis<Shell>(nuclear_framework, basisset_name_beta)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the molecule. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param molecule                     The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name_alpha          The name of the basisset, e.g. "STO-3G", used for the expansion of the alpha spin-orbitals.
     *  @param basisset_name_beta           The name of the basisset, e.g. "STO-3G", used for the expansion of the beta spin-orbitals.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha),
                     ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta)) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create an unrestricted spin-orbital basis from a restricted spin-orbital basis, leading to alpha- and beta- coefficient matrices that are equal.
     *
     *  @param r_spinor_basis               the restricted spinor basis
     *
     *  @return an unrestricted spinor basis
     */
    static USpinorBasis<ExpansionScalar, Shell> FromRestricted(const RSpinorBasis<ExpansionScalar, Shell>& r_spinor_basis) {

        const auto scalar_basis = r_spinor_basis.scalarBasis();
        const auto C = r_spinor_basis.coefficientMatrix();
        return USpinorBasis<ExpansionScalar, Shell>(scalar_basis, scalar_basis, UTransformationMatrix<ExpansionScalar>::FromRestricted(C));
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The transformation that expresses the current spin-orbitals in terms of the underlying scalar basis.
     */
    UTransformationMatrix<ExpansionScalar> coefficientMatrix() const {
        return UTransformationMatrixComponent<ExpansionScalar> {this->alpha().coefficientMatrix(), this->beta().coefficientMatrix()};
    }


    /**
     *  @return The total number of spinors/spin-orbitals that this spin-orbital basis describes.
     */
    size_t numberOfSpinors() const {

        const auto K_alpha = this->alpha().simpleDimension();
        const auto K_beta = this->beta().simpleDimension();

        return K_alpha + K_beta;
    }


    /**
     *  @return The total number of spin-orbitals that this spin-orbital basis describes.
     */
    size_t numberOfSpinOrbitals() const { return this->numberOfSpinors(); }


    /*
     *  MARK: Orthonormality
     */

    /**
     *  Check if this spin-orbital basis is orthonormal, within a given precision.
     * 
     *  @param precision            The precision used to test orthonormality.
     *
     *  @return If this spin-orbital basis is orthonormal.
     */
    bool isOrthonormal(const double precision = 1.0e-08) const { return this->alpha().isOrthonormal(precision) && this->beta().isOrthonormal(precision); }


    /**
     *  @return The transformation T to the Löwdin basis, i.e. T = S_current^{-1/2}.
     */
    UTransformationMatrix<ExpansionScalar> lowdinOrthonormalizationMatrix() const {

        const auto T_a = this->alpha().lowdinOrthonormalizationMatrix();
        const auto T_b = this->beta().lowdinOrthonormalizationMatrix();

        return UTransformationMatrix<ExpansionScalar> {T_a, T_b};
    }


    /**
     *  Transform this spin-orbital basis to the 'Löwdin basis', which is the orthonormal basis characterized by the transformation T = S_current^{-1/2}, where S_current is the current overlap matrix.
     */
    void lowdinOrthonormalize() {
        this->alpha().lowdinOrthonormalize();
        this->beta().lowdinOrthonormalize();
    }


    /*
     *  MARK: Quantizing first-quantized operators
     */


    /**
     *  Quantize a one-electron operator in this unrestricted spin-orbital basis, i.e. express/project the one-electron operator in/onto this spin-orbital basis.
     * 
     *  @param fq_op                                The first-quantized one-electron operator.
     *
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     * 
     *  @return The second-quantized operator corresponding to the given first-quantized operator.
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_op) const -> USQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        // Quantize the one-electron operator in the alpha- and beta- bases and return the wrapped result.
        const auto f_a = this->alpha().quantize(fq_op);
        const auto f_b = this->beta().quantize(fq_op);

        return ResultOperator {f_a, f_b};
    }


    /**
     *  Quantize the Coulomb operator in this unrestricted spin-orbital basis, i.e. express/project the one-electron operator in/onto this spin-orbital basis.
     * 
     *  @param coulomb_op               The first-quantized Coulomb operator operator.
     * 
     *  @return The second-quantized Coulomb operator.
     */
    auto quantize(const CoulombRepulsionOperator& coulomb_op) const -> ScalarUSQTwoElectronOperator<product_t<typename CoulombRepulsionOperator::Scalar, ExpansionScalar>> {

        using ResultScalar = product_t<typename CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = ScalarUSQOneElectronOperator<ResultScalar>;

        // Determine the matrix representation of the four spin-components of the second-quantized Coulomb operator.
        const auto g_aa_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->alpha().scalarBasis(), this->alpha().scalarBasis());
        const auto g_ab_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->alpha().scalarBasis(), this->beta().scalarBasis());
        const auto g_ba_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->beta().scalarBasis(), this->alpha().scalarBasis());
        const auto g_bb_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->beta().scalarBasis(), this->beta().scalarBasis());

        // We have previously calculated the representations in the AO basis, so we'll still have to transform these representations to the current spin-orbitals.
        ResultOperator g {g_aa_par, g_ab_par, g_ba_par, g_bb_par};
        g.transform(this->coefficientMatrix());  // Now, g is expressed in the current spin-orbital basis.

        return g;
    }


    /*
     *  MARK: Mulliken partitioning
     */

    /**
     *  @param ao_list              indices of the AOs used for the atomic spin operator in the z-direction
     *  @param sigma                alpha or beta
     *
     *  @return the SQ atomic spin operator in the z-direction for a set of AOs
     *
     *  @note This method is only available for real SQOperators.
     */
    // template <typename S = ExpansionScalar, typename = IsReal<S>>
    // ScalarUSQOneElectronOperatorComponent<double> calculateAtomicSpinZ(const std::vector<size_t>& ao_list, const Spin& sigma) const {

    //     // The atomic spin operator can be calculated as as the atomic Mulliken operator divided by 2, multiplied by the correct sign factor
    //     int sign = 1 - 2 * sigma;  // 1 for ALPHA, -1 for BETA
    //     const auto spin_z_par = 0.5 * sign * this->spinor_bases[sigma].calculateMullikenOperator(ao_list).parameters();
    //     return ScalarRSQOneElectronOperator<double>(spin_z_par);
    // }


    // /**
    //  *  @param ao_list              indices of the AOs used for the Mulliken populations
    //  *  @param sigma                alpha or beta
    //  *
    //  *  @return the Mulliken operator for a set of AOs and the requested component
    //  *
    //  *  @note This method is only available for real matrix representations.
    //  */
    // template <typename S = ExpansionScalar, typename = IsReal<S>>
    // ScalarRSQOneElectronOperator<double> calculateMullikenOperator(const std::vector<size_t>& ao_list, const Spin& sigma) const { return this->spinor_bases[sigma].template calculateMullikenOperator<ExpansionScalar>(ao_list); }
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<USpinorBasis<_ExpansionScalar, _Shell>> {

    // The type of the transformation matrix for which the basis transformation should be defined. // TODO: Rename "TM" to "TransformationMatrix". A transformation matrix should naturally be transformable with itself.
    using TM = UTransformationMatrix<Scalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<USpinorBasis<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
