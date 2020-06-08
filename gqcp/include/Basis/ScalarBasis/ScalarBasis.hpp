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


#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Mathematical/AbstractFunction/LinearCombination.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"

#include <type_traits>


namespace GQCP {


/**
 *  A class that represents a scalar basis: it represents a collection of scalar basis functions. It provides an interface to obtain basis functions and calculate integrals over the shell type
 *
 * @tparam _Shell       the type of shell that this scalar basis contains
 */
template <typename _Shell>
class ScalarBasis {
public:
    using Shell = _Shell;

    using Primitive = typename Shell::Primitive;          // the type of the primitive functions that underlie this scalar basis
    using BasisFunction = typename Shell::BasisFunction;  // the type of basis functions that this scalar basis consists of


private:
    ShellSet<Shell> shell_set;  // a collection of shells that represents this scalar basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param shell_set        a collection of shells that represents this scalar basis
     */
    ScalarBasis(const ShellSet<Shell>& shell_set) :
        shell_set {shell_set} {}


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework
     *
     *  @param nuclear_framework        the nuclear framework containing the nuclei on which the shells should be centered
     *  @param basisset_name            the name of the basisset, e.g. "STO-3G"
     *
     *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     * 
     *  @note This constructor is only available for GTOShells (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <typename Z = GTOShell>
    ScalarBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(GTOBasisSet(basisset_name).generate(nuclear_framework)) {

        this->shell_set.embedNormalizationFactorsOfPrimitives();
    }


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the molecule
     *
     *  @param molecule             the molecule containing the nuclei on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     *
     *  Note that the normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells
     * 
     *  @note This constructor is only available for GTOShells (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415).
     */
    template <typename Z = GTOShell>
    ScalarBasis(const Molecule& molecule, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(molecule.nuclearFramework(), basisset_name) {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the basis functions that 'are' in this scalar basis
     */
    std::vector<BasisFunction> basisFunctions() const { return this->shell_set.basisFunctions(); }

    /**
     *  @return the number of basis functions that 'are' in this scalar basis
     */
    size_t numberOfBasisFunctions() const { return this->shell_set.numberOfBasisFunctions(); }

    /**
     *  @return the underlying set of shells
     */
    const ShellSet<Shell>& shellSet() const { return this->shell_set; }
};


}  // namespace GQCP
