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
#ifndef ShellSet_hpp
#define ShellSet_hpp


#include "Basis/Shell.hpp"
#include "Molecule.hpp"

#include <vector>


namespace GQCP {


/**
 *  A class that represents a list of shells (and therefore extends std::vector<Shell>)
 */
class ShellSet : public std::vector<Shell> {
public:
    using std::vector<Shell>::vector;  // inherit base constructors


public:
    // CONSTRUCTORS
    ShellSet() = default;
    /**
     *  Construct a ShellSet by placing the shells corresponding to the basisset information on every atom of the molecule
     *
     *  @param molecule             the molecule containing the atoms on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     */
    ShellSet(const Molecule& molecule, const std::string& basisset_name);


    // PUBLIC METHODS
    /**
     *  @return the number of shells in this shell set
     */
    size_t numberOfShells() const;

    /**
     *  @return the number of basis functions in this shell set
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  @return an ordered vector of the unique atoms in this shell set
     */
    std::vector<Atom> atoms() const;

    /**
     *  @param shell_index      the index of the shell
     *
     *  @return the (total basis function) index that corresponds to the first basis function in the given shell
     */
    size_t basisFunctionIndex(size_t shell_index) const;

    /**
     *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void embedNormalizationFactorsOfPrimitives();

    /**
     *  For every of the shells, embed the normalization factor of every Gaussian primitive into its corresponding contraction coefficient. If this has already been done, this function does nothing
     *
     *  Note that the normalization factor that is embedded corresponds to the spherical (or axis-aligned Cartesian) GTO
     */
    void unEmbedNormalizationFactorsOfPrimitives();

    /**
     *  For every of the shells, embed the total normalization factor of the corresponding linear combination of spherical (or axis-aligned Cartesian) GTOs into the contraction coefficients
     */
    void embedNormalizationFactors();
};


}  // namespace GQCP


#endif  /* ShellSet_hpp */
