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
#include "Basis/ShellSet.hpp"

#include "Basis/LibintInterfacer.hpp"

#include <algorithm>


namespace GQCP {



/*
 *  CONSTRUCTORS
 */

/**
 *  Construct a ShellSet by placing the shells corresponding to the basisset information on every nucleus of the molecule
 *
 *  @param molecule             the molecule containing the nuclei on which the shells should be centered
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 */
ShellSet::ShellSet(const Molecule& molecule, const std::string& basisset_name) {

    // Since we haven't implemented our own BasisSet class, we use libint to read in the basisset specification file and create the shells
    // TODO no longer use libint2 to read this

    const auto& nuclei = molecule.nuclearFramework().nucleiAsVector();

    libint2::BasisSet libint_basis (basisset_name, LibintInterfacer::get().interface(nuclei));

    *this = LibintInterfacer::get().interface(libint_basis, nuclei);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of shells in this shell set
 */
size_t ShellSet::numberOfShells() const {
    return this->size();
}


/**
 *  @return the number of basis functions in this shell set
 */
size_t ShellSet::numberOfBasisFunctions() const {

    size_t nbf {};
    for (const auto& shell : *this) {
        nbf += shell.numberOfBasisFunctions();
    }
    return nbf;
}


/**
 *  @return an ordered vector of the unique nuclei in this shell set
 */
std::vector<Nucleus> ShellSet::nuclei() const {

    std::vector<Nucleus> nuclei {};

    // Append every unique nucleus in this shell set's shells
    for (const auto& shell : *this) {
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
size_t ShellSet::basisFunctionIndex(size_t shell_index) const {

    size_t bf_index {};

    // Count the number of basis functions before
    for (size_t i = 0; i < shell_index; i++) {
        bf_index += this->operator[](i).numberOfBasisFunctions();
    }

    return bf_index;
}


/**
 *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
 *
 *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
 */
void ShellSet::embedNormalizationFactorsOfPrimitives() {

    for (auto& shell : *this) {
        shell.embedNormalizationFactorsOfPrimitives();
    }
}

    
/**
 *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
 *
 *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
 */
void ShellSet::unEmbedNormalizationFactorsOfPrimitives() {

    for (auto& shell : *this) {
        shell.unEmbedNormalizationFactorsOfPrimitives();
    }
}


/**
 *  For every of the shells, embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients
 */
void ShellSet::embedNormalizationFactors() {

    for (auto& shell : *this) {
        shell.embedNormalizationFactor();
    }
}


}  // namespace GQCP
