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


#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"


namespace GQCP {


/**
 *  A type specifically designed to act as a parent class for `RSpinOrbitalBasis` and `USpinOrbitalBasisComponent` in order to share common functionality.
 * 
 *  @tparam _ExpansionScalar                The scalar type used to represent an expansion coefficient of the spin-orbitals in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                          The type of shell that the underlying scalar basis contains.
 *  @tparam _DerivedSpinOrbitalBasis        The spin-orbital basis that derives from this class, enabling CRTP and compile-time polymorphism.
 */
template <typename _ExpansionScalar, typename _Shell, typename _DerivedSpinOrbitalBasis>
class SimpleSpinOrbitalBasis:
    public SimpleSpinorBasis<_ExpansionScalar, _DerivedSpinOrbitalBasis> {
public:
    // The scalar type used to represent an expansion coefficient of the spin-orbitals in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell that the underlying scalar basis contains.
    using Shell = _Shell;

    // The spin-orbital basis that derives from this class, enabling CRTP and compile-time polymorphism.
    using DerivedSpinOrbitalBasis = _DerivedSpinOrbitalBasis;

    // The type of the base spinor basis.
    using BaseSpinorBasis = SimpleSpinorBasis<ExpansionScalar, DerivedSpinOrbitalBasis>;

    // The type of transformation that is naturally related to the derived spin-orbital basis.
    using Transformation = typename SpinorBasisTraits<DerivedSpinOrbitalBasis>::Transformation;


protected:
    // The underlying scalar basis with respect to which the basis coefficients are expressed.
    ScalarBasis<Shell> scalar_basis;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a `SimpleSpinOrbitalBasis` from a scalar basis and a transformation that expresses the current spin-orbitals in terms of that underlying scalar basis.
     * 
     *  @param scalar_basis         The underlying scalar basis with respect to which the basis coefficients are expressed.
     *  @param C                    The transformation that relates the current set of spinors with the atomic spinors.
     */
    SimpleSpinOrbitalBasis(const ScalarBasis<Shell>& scalar_basis, const Transformation& C) :
        BaseSpinorBasis(C),
        scalar_basis {scalar_basis} {}


    /**
     *  Construct a `SimpleSpinOrbitalBasis` with an initial coefficient matrix that is the identity. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     * 
     *  @param scalar_basis         The underlying scalar basis with respect to which the basis coefficients are expressed.
     */
    SimpleSpinOrbitalBasis(const ScalarBasis<Shell>& scalar_basis) :
        SimpleSpinOrbitalBasis(scalar_basis, Transformation::Identity(scalar_basis.numberOfBasisFunctions())) {}


    /**
     *  Construct a `SimpleSpinOrbitalBasis` with an underlying scalar basis that is made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells of the scalar basis should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    SimpleSpinOrbitalBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        SimpleSpinOrbitalBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct a `SimpleSpinOrbitalBasis` with an underlying scalar basis that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spin-orbital basis then corresponds to the (non-orthogonal) atomic spin-orbitals (AOs).
     *
     *  @param molecule             The molecule containing the nuclei on which the shells of the scalar basis should be centered.
     *  @param basisset_name        The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    SimpleSpinOrbitalBasis(const Molecule& molecule, const std::string& basisset_name) :
        SimpleSpinOrbitalBasis(molecule.nuclearFramework(), basisset_name) {}


    /**
     *  Construct a simple spin-orbital basis with an underlying scalar basis that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spinor basis corresponds to the non-orthogonal London atomic spinors (AOs).
     *
     *  @param molecule                 The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *  @param B                        The homogeneous magnetic field.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    SimpleSpinOrbitalBasis(const Molecule& molecule, const std::string& basisset_name, const HomogeneousMagneticField& B,
                           typename std::enable_if<std::is_same<Z, LondonGTOShell>::value>::type* = 0) :
        SimpleSpinOrbitalBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name, B)) {}


    /*
     *  MARK: Scalar basis
     */

    /**
     *  @return The underlying scalar basis with respect to which the basis coefficients are expressed.
     */
    const ScalarBasis<Shell>& scalarBasis() const { return this->scalar_basis; }
};


}  // namespace GQCP
