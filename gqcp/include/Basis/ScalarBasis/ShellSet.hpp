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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Molecule/Molecule.hpp"
#include "Physical/HomogeneousMagneticField.hpp"
#include "Utilities/type_traits.hpp"

#include <algorithm>
#include <initializer_list>
#include <vector>


namespace GQCP {


/**
 *  A collection of shells.
 * 
 *  @tparam _Shell          The type of shell that this shell set contains.
 */
template <typename _Shell>
class ShellSet {
public:
    // The type of shell that this shell set contains.
    using Shell = _Shell;

    // The type of primitives that this shell set contains.
    using Primitive = typename Shell::Primitive;

    // The type of basis functions that this shell contains.
    using BasisFunction = typename Shell::BasisFunction;


private:
    // The collection of shells represented by a vector.
    std::vector<Shell> shells;


public:
    /*
     *  MARK: Constructors
     */

    /**
     * @param shells            The collection of shells represented by a vector.
     */
    ShellSet(const std::vector<Shell>& shells) :
        shells {shells} {}


    /**
     *  Construct a ShellSet using an `initializer_list`.
     * 
     *  @param list             The `initializer_list`.
     */
    ShellSet(const std::initializer_list<GTOShell>& list) :
        shells {list} {}


    /*
     *  MARK: Normalization
     */

    /**
     *  For every of the shells, embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients.
     */
    void embedNormalizationFactors() {

        for (auto& shell : this->shells) {
            shell.embedNormalizationFactor();
        }
    }


    /**
     *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function has no effect.
     *
     *  @note The normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO.
     */
    void embedNormalizationFactorsOfPrimitives() {

        for (auto& shell : this->shells) {
            shell.embedNormalizationFactorsOfPrimitives();
        }
    }


    /**
     *  For every of the shells, unembed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function has no effect.
     *
     *  @note The normalization factor that is unembedded corresponds to the spherical (or axis-aligned Cartesian) GTO.
     */
    void unEmbedNormalizationFactorsOfPrimitives() {

        for (auto& shell : this->shells) {
            shell.unEmbedNormalizationFactorsOfPrimitives();
        }
    }


    /*
     *  MARK: Basis functions
     */

    /**
     *  @return The number of basis functions in this shell set.
     */
    size_t numberOfBasisFunctions() const {

        size_t value {};
        for (const auto& shell : this->shells) {
            value += shell.numberOfBasisFunctions();
        }

        return value;
    }


    /**
     *  @return The basis functions that are in this shell.
     */
    std::vector<BasisFunction> basisFunctions() const {

        std::vector<BasisFunction> basis_functions;
        basis_functions.reserve(this->numberOfBasisFunctions());
        for (const auto& shell : this->shells) {
            const auto shell_basis_functions = shell.basisFunctions();
            for (const auto& shell_basis_function : shell_basis_functions) {
                basis_functions.push_back(shell_basis_function);
            }
        }

        return basis_functions;
    }


    /**
     *  @param shell_index      The index of the shell.
     *
     *  @return The index of the given shell's first basis function in the total set of basis functions of this shell.
     */
    size_t basisFunctionIndex(const size_t shell_index) const {

        // Count the number of basis functions before the given index.
        size_t bf_index {};
        for (size_t i = 0; i < shell_index; i++) {
            bf_index += this->shells[i].numberOfBasisFunctions();
        }

        return bf_index;
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return The collection of this shell set's shells represented by a vector.
     */
    const std::vector<Shell>& asVector() const { return this->shells; }

    /**
     *  @return The number of shells in this shell set.
     */
    size_t numberOfShells() const { return this->shells.size(); }

    /**
     *  @return An ordered vector of the unique nuclei in this shell set.
     */
    std::vector<Nucleus> nuclei() const {

        // Append every unique nucleus in this shell set's shells.
        std::vector<Nucleus> nuclei {};
        for (const auto& shell : this->shells) {
            const auto& nucleus = shell.nucleus();

            const auto unary_predicate = [nucleus](const Nucleus& other) {
                return Nucleus::equalityComparer()(nucleus, other);
            };
            const auto& p = std::find_if(nuclei.begin(), nuclei.end(), unary_predicate);

            if (p == nuclei.end()) {  // If the nucleus is unique.
                nuclei.push_back(nucleus);
            }
        }

        return nuclei;
    }


    /**
     *  @return The maximum angular momentum of the shells.
     */
    size_t maximumAngularMomentum() const {

        // Generate a list of all the angular momenta and then check its maximum.
        std::vector<size_t> angular_momenta {};  // This will contain the number of primitives for each of the shells.
        angular_momenta.reserve(this->numberOfShells());

        for (const auto& shell : this->asVector()) {
            angular_momenta.push_back(shell.angularMomentum());
        }

        const auto it = std::max_element(angular_momenta.begin(), angular_momenta.end());  // 'it' for iterator.
        return *it;
    }


    /**
     *  @return The maximum number of primitives that are used inside the shells.
     */
    size_t maximumNumberOfPrimitives() const {

        // Generate a list of all the number of primitives momenta and then check its maximum.
        std::vector<size_t> number_of_primitives {};  // This will contain the number of primitives for each of the shells.
        number_of_primitives.reserve(this->numberOfShells());

        for (const auto& shell : this->asVector()) {
            number_of_primitives.push_back(shell.contractionSize());
        }

        const auto it = std::max_element(number_of_primitives.begin(), number_of_primitives.end());  // 'it' for iterator.
        return *it;
    }


    /*
     *  MARK: London modifications
     */

    /**
     *  Apply a London- (gauge-including) modification to each shell.
     * 
     *  @param B            The homogeneous magnetic field.
     * 
     *  @return A new shell shet whose underlying shells have been London-modified.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, GTOShell>::value, ShellSet<LondonGTOShell>> applyLondonModification(const HomogeneousMagneticField& B) const {

        // Run over each `GTOShell` and construct the corresponding `LondonGTOShell`.
        std::vector<LondonGTOShell> london_shells;
        london_shells.reserve(this->numberOfShells());
        for (const auto& gto_shell : this->asVector()) {
            LondonGTOShell london_shell {gto_shell, B};
            london_shells.push_back(london_shell);
        }

        return ShellSet<LondonGTOShell>(london_shells);
    }
};


}  // namespace GQCP
