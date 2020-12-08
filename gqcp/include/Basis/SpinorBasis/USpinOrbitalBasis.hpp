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


#include "Basis/MullikenPartitioning/UMullikenPartitioning.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasisComponent.hpp"
#include "Basis/Transformations/SpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/SpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
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
class USpinOrbitalBasis:
    public SpinResolvedBase<USpinOrbitalBasisComponent<_ExpansionScalar, _Shell>, USpinOrbitalBasis<_ExpansionScalar, _Shell>>,
    public SpinResolvedBasisTransformable<USpinOrbitalBasis<_ExpansionScalar, _Shell>>,
    public SpinResolvedJacobiRotatable<USpinOrbitalBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type of 'this'.
    using Self = USpinOrbitalBasis<ExpansionScalar, Shell>;

    // The type that is used for representing the primitive for a basis function of this spin-orbital basis' underlying AO basis.
    using Primitive = typename Shell::Primitive;

    // The type that is used for representing the underlying basis functions of this spin-orbital basis.
    using BasisFunction = typename Shell::BasisFunction;

    // The type component this spin resolved object is made of.
    using ComponentType = typename SpinResolvedBase<USpinOrbitalBasisComponent<ExpansionScalar, Shell>, Self>::Of;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<USpinOrbitalBasisComponent<ExpansionScalar, Shell>, USpinOrbitalBasis<ExpansionScalar, Shell>>::SpinResolvedBase;


    /**
     *  Create an `USpinOrbitalBasis` from an alpha and a beta-spin-orbital basis and a transformation that expresses the current spin-orbitals in terms of that underlying scalar basis.
     * 
     *  @param alpha_scalar_basis           The scalar basis in which the alpha spin-orbitals are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta spin-orbitals are expanded.
     *  @param C                            The transformation that relates the current set of spinors with the atomic spinors.
     */
    USpinOrbitalBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const UTransformation<ExpansionScalar>& C) :
        USpinOrbitalBasis(USpinOrbitalBasisComponent<ExpansionScalar, Shell> {alpha_scalar_basis, C.alpha()},
                          USpinOrbitalBasisComponent<ExpansionScalar, Shell> {beta_scalar_basis, C.beta()}) {

        // Check if the dimensions of the given objects are compatible.
        if (C.alpha().numberOfOrbitals() != alpha_scalar_basis.numberOfBasisFunctions()) {
            throw std::invalid_argument("USpinOrbitalBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const UTransformation<ExpansionScalar>&): The given dimensions of the scalar basis and coefficient matrix for the alpha spin-orbitals are incompatible.");
        }

        if (C.beta().numberOfOrbitals() != beta_scalar_basis.numberOfBasisFunctions()) {
            throw std::invalid_argument("USpinOrbitalBasis(const ScalarBasis<Shell>&, const ScalarBasis<Shell>&, const UTransformation<ExpansionScalar>&): The given dimensions of the scalar basis and coefficient matrix for the beta spin-orbitals are incompatible.");
        }
    }


    /**
     *  Construct an unrestricted spin-orbital basis in which both underlying scalar bases and their expansions are equal.
     *
     *  @param scalar_basis         The scalar basis in which both the alpha and beta spin-orbitals are expanded.
     *  @param C                    The transformation that expresses the current spin-orbitals in terms of the underlying scalar basis.
     */
    USpinOrbitalBasis(const ScalarBasis<Shell>& scalar_basis, const UTransformationComponent<ExpansionScalar>& C) :
        USpinOrbitalBasis(scalar_basis, scalar_basis, UTransformation<ExpansionScalar>::FromEqual(C)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with two different underlying scalar bases, and a coefficient matrix being the identity. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     * 
     *  @param alpha_scalar_basis           The scalar basis in which the alpha spin-orbitals are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta spin-orbitals are expanded.
     */
    USpinOrbitalBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        USpinOrbitalBasis(alpha_scalar_basis, beta_scalar_basis,
                          UTransformation<ExpansionScalar>(
                              UTransformationComponent<ExpansionScalar>::Identity(alpha_scalar_basis.numberOfBasisFunctions()),
                              UTransformationComponent<ExpansionScalar>::Identity(beta_scalar_basis.numberOfBasisFunctions()))) {}


    /**
     *  Construct an unrestricted spin-orbital basis in which both underlying scalar bases are equal, and the coefficient matrix being the identity. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param scalar_basis         The scalar basis in which both the alpha and beta spin-orbitals are expanded.
     */
    USpinOrbitalBasis(const ScalarBasis<Shell>& scalar_basis) :
        USpinOrbitalBasis(scalar_basis, scalar_basis) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases (equal for both the alpha and beta components) that are made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells of the scalar basis should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinOrbitalBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        USpinOrbitalBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases (equal for both the alpha and beta components) that are made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param molecule                 The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinOrbitalBasis(const Molecule& molecule, const std::string& basisset_name) :
        USpinOrbitalBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name)) {}


    /**
     *  Construct an unrestricted spin-orbital basis with underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the nuclear framework. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param nuclear_framework            The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name_alpha          The name of the basisset, e.g. "STO-3G", used for the expansion of the alpha spin-orbitals.
     *  @param basisset_name_beta           The name of the basisset, e.g. "STO-3G", used for the expansion of the beta spin-orbitals.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    USpinOrbitalBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinOrbitalBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha),
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
    USpinOrbitalBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        USpinOrbitalBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha),
                          ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta)) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Create an unrestricted spin-orbital basis from a restricted spin-orbital basis, leading to alpha- and beta- coefficient matrices that are equal.
     *
     *  @param r_spinor_basis               The restricted spin-orbital basis.
     *
     *  @return An `USpinOrbitalBasis` that corresponds to the given restricted one.
     */
    static USpinOrbitalBasis<ExpansionScalar, Shell> FromRestricted(const RSpinOrbitalBasis<ExpansionScalar, Shell>& r_spinor_basis) {

        const auto& scalar_basis = r_spinor_basis.scalarBasis();
        const auto& C = r_spinor_basis.expansion();

        return USpinOrbitalBasis<ExpansionScalar, Shell>(scalar_basis, scalar_basis, UTransformation<ExpansionScalar>::FromRestricted(C));
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The transformation that expresses the current spin-orbitals in terms of the underlying scalar basis.
     */
    UTransformation<ExpansionScalar> expansion() const {
        return UTransformation<ExpansionScalar> {this->alpha().expansion(), this->beta().expansion()};
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
     *  @return The transformation to the Löwdin basis. See also `SimpleSpinOrbitalBasis`.
     */
    UTransformation<ExpansionScalar> lowdinOrthonormalization() const {

        const auto T_a = this->alpha().lowdinOrthonormalization();
        const auto T_b = this->beta().lowdinOrthonormalization();

        return UTransformation<ExpansionScalar> {T_a, T_b};
    }


    /**
     *  Transform this spin-orbital basis to the 'Löwdin basis'. See also `SimpleSpinOrbitalBasis`.
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
     *  Quantize the z-component of the electronic spin operator in this unrestricted spin-orbital basis, i.e. express/project the one-electron operator in/onto this spin-orbital basis.
     * 
     *  @param fq_op                                The first-quantized one-electron operator.
     *
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     * 
     *  @return The second-quantized operator corresponding to the given first-quantized operator.
     */
    auto quantize(const ElectronicSpin_zOperator& fq_op) const -> USQOneElectronOperator<product_t<ElectronicSpin_zOperator::Scalar, ExpansionScalar>, ElectronicSpin_zOperator::Vectorizer> {

        using ResultScalar = product_t<ElectronicSpin_zOperator::Scalar, ExpansionScalar>;
        using ResultOperator = USQOneElectronOperator<ResultScalar, ElectronicSpin_zOperator::Vectorizer>;

        // We can use the quantization of the overlap operator for the quantization of S_z.
        const auto S_a = this->alpha().overlap();  // In the current orbital basis.
        const auto S_b = this->beta().overlap();   // In the current orbital basis.

        const auto S_z_a = 0.5 * S_a;
        const auto S_z_b = -0.5 * S_b;

        return ResultOperator {S_z_a, S_z_b};
    }


    /**
     *  Quantize the overlap operator in this spin-orbital basis.
     * 
     *  @return The second-quantized overlap operator.
     */
    ScalarUSQOneElectronOperator<ExpansionScalar> overlap() const { return this->quantize(Operator::Overlap()); }


    /**
     *  Quantize the Coulomb operator in this unrestricted spin-orbital basis, i.e. express/project the one-electron operator in/onto this spin-orbital basis.
     * 
     *  @param coulomb_op               The first-quantized Coulomb operator operator.
     * 
     *  @return The second-quantized Coulomb operator.
     */
    auto quantize(const CoulombRepulsionOperator& coulomb_op) const -> ScalarUSQTwoElectronOperator<product_t<typename CoulombRepulsionOperator::Scalar, ExpansionScalar>> {

        using ResultScalar = product_t<typename CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = ScalarUSQTwoElectronOperator<ResultScalar>;

        // Determine the matrix representation of the four spin-components of the second-quantized Coulomb operator.
        const auto g_aa_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->alpha().scalarBasis(), this->alpha().scalarBasis());
        const auto g_ab_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->alpha().scalarBasis(), this->beta().scalarBasis());
        const auto g_ba_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->beta().scalarBasis(), this->alpha().scalarBasis());
        const auto g_bb_par = IntegralCalculator::calculateLibintIntegrals(coulomb_op, this->beta().scalarBasis(), this->beta().scalarBasis());

        // We have previously calculated the representations in the AO basis, so we'll still have to transform these representations to the current spin-orbitals.
        ResultOperator g {g_aa_par, g_ab_par, g_ba_par, g_bb_par};
        g.transform(this->expansion());  // Now, g is expressed in the current spin-orbital basis.

        return g;
    }


    /**
     *  MARK: Mulliken partitioning
     */

    /**
     *  Partition this set of unrestricted spin-orbitals according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken partitioning.
     * 
     *  @return A `UMullikenPartitioning` for the AOs selected by the supplied selector function.
     */
    UMullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const BasisFunction&)>& selector) const {

        return UMullikenPartitioning<ExpansionScalar> {this->alpha().mullikenPartitioning(selector), this->beta().mullikenPartitioning(selector)};
    }


    /**
     *  Partition this set of unrestricted spin-orbitals according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for shells that should be included the Mulliken partitioning.
     * 
     *  @return A `UMullikenPartitioning` for the AOs selected by the supplied selector function.
     */
    UMullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const Shell&)>& selector) const {

        return UMullikenPartitioning<ExpansionScalar> {this->alpha().mullikenPartitioning(selector), this->beta().mullikenPartitioning(selector)};
    }


    /**
     *  MARK: Enabling basis transformations
     */

    // Since `rotate` and `rotated` are both defined in `SpinResolvedBasisTransformable` and `SpinResolvedJacobiRotatable`, we have to explicitly enable these methods here.

    // Allow the `rotate` method from `SpinResolvedBasisTransformable`, since there's also a `rotate` from `SpinResolvedJacobiRotatable`.
    using SpinResolvedBasisTransformable<Self>::rotate;

    // Allow the `rotated` method from `SpinResolvedBasisTransformable`, since there's also a `rotated` from `SpinResolvedJacobiRotatable`.
    using SpinResolvedBasisTransformable<Self>::rotated;

    // Allow the `rotate` method from `SpinResolvedJacobiRotatable`, since there's also a `rotate` from `SpinResolvedBasisTransformable`.
    using SpinResolvedJacobiRotatable<Self>::rotate;

    // Allow the `rotated` method from `SpinResolvedJacobiRotatable`, since there's also a `rotated` from `SpinResolvedBasisTransformable`.
    using SpinResolvedJacobiRotatable<Self>::rotated;
};

/*
 *  MARK: Convenience aliases
 */
template <typename ExpansionScalar, typename Shell>
using USpinorBasis = USpinOrbitalBasis<ExpansionScalar, Shell>;


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<USpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of transformation that is naturally related to a `USpinOrbitalBasis`.
    using Transformation = UTransformation<_ExpansionScalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<USpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation for which the Jacobi rotation should be defined.
    using JacobiRotationType = UJacobiRotation;
};


}  // namespace GQCP
