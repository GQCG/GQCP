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


#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Molecule/Molecule.hpp"

#include <algorithm>
#include <initializer_list>
#include <vector>


namespace GQCP {


/**
 *  A class that represents a set of shells
 * 
 *  @tparam _Shell          the type of shell that is contained in this set
 */
template <typename _Shell>
class ShellSet {
public:
    using Shell = _Shell;  // the type of shell that is contained in this set
    using BasisFunction = typename Shell::BasisFunction;  // the type of basis function that the shell can produce


private:
    std::vector<Shell> shells;  // all the shells represented by a vector


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     * @param shells            all the shells represented by a vector
     */
    ShellSet(const std::vector<Shell>& shells) : 
        shells (shells)
    {}


    /**
     *  Construct a ShellSet using an initializer_list
     */
    ShellSet(const std::initializer_list<GTOShell>& list) :
        shells (list)
    {}



    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return a vector of the underlying shells that this shellset describes
     */
    const std::vector<Shell>& asVector() const { return this->shells; }

    /**
     *  @return the maximum angular momentum of the shells
     */
    size_t maximumAngularMomentum() const {

        std::vector<size_t> angular_momenta (this->numberOfShells());  // contains the number of primitives for each of the shells

        for (size_t i = 0; i < this->numberOfShells(); i++) {
            const auto shell = this->shells[i];

            angular_momenta[i] = shell.angularMomentum();
        }

        const auto it = std::max_element(angular_momenta.begin(), angular_momenta.end());  // iterator
        return *it;
    }


    /**
     *  @return the maximum number of primitives that are used inside the shells
     */
    size_t maximumNumberOfPrimitives() const {

        std::vector<size_t> number_of_primitives (this->numberOfShells());  // contains the number of primitives for each of the shells

        for (size_t i = 0; i < this->numberOfShells(); i++) {
            const auto shell = this->shells[i];

            number_of_primitives[i] = shell.contractionSize();
        }

        const auto it = std::max_element(number_of_primitives.begin(), number_of_primitives.end());  // iterator
        return *it;
    }
    

    /**
     *  @return the number of shells in this shell set
     */
    size_t numberOfShells() const { return this->shells.size(); }

    /**
     *  @return the number of basis functions in this shell set
     */
    size_t numberOfBasisFunctions() const {

        size_t value {};
        for (const auto& shell : this->shells) {
            value += shell.numberOfBasisFunctions();
        }

        return value;
    }


    /**
     *  @return an ordered vector of the unique nuclei in this shell set
     */
    std::vector<Nucleus> nuclei() const {

        // Append every unique nucleus in this shell set's shells
        std::vector<Nucleus> nuclei {};
        for (const auto& shell : this->shells) {
            const auto& nucleus = shell.get_nucleus();

            const auto unary_predicate = [nucleus] (const Nucleus& other) {
                return Nucleus::equalityComparer()(nucleus, other);
            };
            const auto& p = std::find_if(nuclei.begin(), nuclei.end(), unary_predicate);

            if (p == nuclei.end()) {  // if unique
                nuclei.push_back(nucleus);
            }
        }

        return nuclei;
    }


    /**
     *  @param shell_index      the index of the shell
     *
     *  @return the (total basis function) index that corresponds to the first basis function in the given shell
     */
    size_t basisFunctionIndex(size_t shell_index) const {

        // Count the number of basis functions before the given index
        size_t bf_index {};
        for (size_t i = 0; i < shell_index; i++) {
            bf_index += this->shells[i].numberOfBasisFunctions();
        }

        return bf_index;
    }


    /**
     *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void embedNormalizationFactorsOfPrimitives() {

        for (auto& shell : this->shells) {
            shell.embedNormalizationFactorsOfPrimitives();
        }
    }


    /**
     *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void unEmbedNormalizationFactorsOfPrimitives() {

        for (auto& shell : this->shells) {
            shell.unEmbedNormalizationFactorsOfPrimitives();
        }
    }


    /**
     *  For every of the shells, embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients
     */
    void embedNormalizationFactors() {

        for (auto& shell : this->shells) {
            shell.embedNormalizationFactor();
        }
    }

    /**
     *  @return the basis functions that 'are' in this shell
     */
    std::vector<LinearCombination<double, BasisFunction>> basisFunctions() const {

        throw std::runtime_error("ShellSet::basisFunctions(): This method has not been implemented yet.");
    }
};


}  // namespace GQCP
