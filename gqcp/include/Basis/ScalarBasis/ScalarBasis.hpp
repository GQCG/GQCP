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
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Utilities/type_traits.hpp"

#include <functional>


namespace GQCP {


/**
 *  A class that represents a scalar basis: it represents a collection of scalar basis functions. It provides an interface to obtain basis functions and calculate integrals over the shell type.
 *
 * @tparam _Shell           The type of shell that this scalar basis contains.
 */
template <typename _Shell>
class ScalarBasis {
public:
    // The type of shell that this scalar basis contains.
    using Shell = _Shell;

    // The type of the primitive functions that underlie this scalar basis.
    using Primitive = typename Shell::Primitive;

    // The type of basis functions that this scalar basis consists of.
    using BasisFunction = typename Shell::BasisFunction;

private:
    // A collection of shells that represents this scalar basis.
    ShellSet<Shell> shell_set;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param shell_set        A collection of shells that represents this scalar basis.
     */
    ScalarBasis(const ShellSet<Shell>& shell_set) :
        shell_set {shell_set} {}


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework.
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    ScalarBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(GTOBasisSet(basisset_name).generate(nuclear_framework)) {

        this->shell_set.embedNormalizationFactorsOfPrimitives();
    }


    /**
     *  Construct a scalar basis by placing shells corresponding to the basisset specification on every nucleus of the molecule.
     *
     *  @param molecule             The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name        The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    ScalarBasis(const Molecule& molecule, const std::string& basisset_name,
                typename std::enable_if<std::is_same<Z, GTOShell>::value>::type* = 0) :
        ScalarBasis(molecule.nuclearFramework(), basisset_name) {}


    /**
     *  Construct a scalar basis by placing gauge-including shells corresponding to the basisset specification on every nucleus of the nuclear framework.
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *  @param B                        The homogeneous magnetic field.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    ScalarBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name, const HomogeneousMagneticField& B,
                typename std::enable_if<std::is_same<Z, LondonGTOShell>::value>::type* = 0) :
        ScalarBasis(GTOBasisSet(basisset_name).generate(nuclear_framework).applyLondonModification(B)) {

        this->shell_set.embedNormalizationFactorsOfPrimitives();
    }


    /**
     *  Construct a scalar basis by placing gauge-including shells corresponding to the basisset specification on every nucleus of the molecule.
     *
     *  @param molecule             The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name        The name of the basisset, e.g. "STO-3G".
     *  @param B                    The homogeneous magnetic field.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    ScalarBasis(const Molecule& molecule, const std::string& basisset_name, const HomogeneousMagneticField& B,
                typename std::enable_if<std::is_same<Z, LondonGTOShell>::value>::type* = 0) :
        ScalarBasis(molecule.nuclearFramework(), basisset_name, B) {}


    /*
     *  MARK: Shell set
     */

    /**
     *  @return The collection of shells that represents this scalar basis.
     */
    const ShellSet<Shell>& shellSet() const { return this->shell_set; }


    /*
     *  MARK: Basis functions
     */

    /**
     *  @return The basis functions that this scalar basis consists of.
     */
    std::vector<BasisFunction> basisFunctions() const { return this->shell_set.basisFunctions(); }

    /**
     *  @return The number of basis functions that this scalar basis consists of.
     */
    size_t numberOfBasisFunctions() const { return this->shell_set.numberOfBasisFunctions(); }

    /**
     *  Find the basis function indices selected by a given selector.
     * 
     *  @param selector             A function that returns true if the basis function indices of a particular shell should be included.
     * 
     *  @return The indices of the basis functions whose corresponding shells are selected.
     */
    std::vector<size_t> basisFunctionIndices(const std::function<bool(const Shell&)>& selector) const {

        const auto shells = this->shellSet().asVector();

        // Find the indices of those basis functions for which the shell selector returns true.
        std::vector<size_t> ao_indices;
        size_t bf_index = 0;
        for (size_t shell_index = 0; shell_index < shells.size(); shell_index++) {
            const auto& shell = shells[shell_index];

            // If a shell has to be included, include all indices of the basis functions in it.
            const auto number_of_bf_in_shell = shell.numberOfBasisFunctions();
            if (selector(shell)) {
                for (size_t i = 0; i < number_of_bf_in_shell; i++) {
                    ao_indices.push_back(bf_index);
                    bf_index++;
                }
            } else {
                // Increase the current BF index to accommodate to the next shell.
                bf_index += number_of_bf_in_shell;
            }
        }

        return ao_indices;
    }


    /**
     *  Find the basis function indices selected by a given selector.
     * 
     *  @param selector             A function that returns true if the basis function index should be included.
     * 
     *  @return The indices of the basis functions whose corresponding shells are selected.
     */
    std::vector<size_t> basisFunctionIndices(const std::function<bool(const BasisFunction&)>& selector) const {

        const auto basis_functions = this->basisFunctions();

        // Find the indices of those basis functions for which the selector returns true.
        std::vector<size_t> ao_indices;
        for (size_t i = 0; i < basis_functions.size(); i++) {
            const auto& basis_function = basis_functions[i];

            if (selector(basis_function)) {
                ao_indices.push_back(i);
            }
        }

        return ao_indices;
    }
};


}  // namespace GQCP
