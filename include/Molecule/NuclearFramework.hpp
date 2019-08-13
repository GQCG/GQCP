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
#ifndef GQCP_NUCLEARFRAMEWORK_HPP
#define GQCP_NUCLEARFRAMEWORK_HPP


#include "Mathematical/CartesianDirection.hpp"
#include "Molecule/Nucleus.hpp"

#include <vector>



namespace GQCP {


/**
 *  A class that represents a set of nuclei at fixed positions (in bohr) in space
 */
class NuclearFramework {
private:
    std::vector<Nucleus> nuclei;  // coordinates in bohr


public:
    // CONSTRUCTORS

    /**
     *  @param nuclei           the nuclei of the nuclear framework
     */
    NuclearFramework(const std::vector<Nucleus>& nuclei);


    // NAMED CONSTRUCTORS

    /**
     *  Construct a nuclear framework based on the content of a given .xyz-file. In an .xyz-file, the nuclear coordinates are in Angstrom
     *
     *  @param xyz_filename     the .xyz-file that contains the nuclear coordinates in Angstrom
     */
    static NuclearFramework ReadXYZ(const std::string& xyz_filename);

    /**
     *  @param n            the number of H nuclei
     *  @param spacing      the internuclear spacing in bohr
     *  @param axis         the Cartesian axis on which the H-chain should be placed
     *
     *  @return a H-chain with equal internuclear spacing
     */
    static NuclearFramework HChain(const size_t n, const double spacing, CartesianDirection axis=CartesianDirection::z);

    /**
     *  @param n        the number of H2-molecules
     *  @param a        the internuclear distance in bohr
     *  @param b        the intermolecular distance in bohr
     *  @param axis     the Cartesian axis on which the H2-chain should be placed
     * 
     *  @return a H2-chain with the specified internuclear and intermolecular distances
     */
    static NuclearFramework H2Chain(const size_t n, const double a, const double b, CartesianDirection axis=CartesianDirection::z);


    // OPERATORS

    /**
     *  @param os                       the output stream which the nuclear framework should be concatenated to
     *  @param nuclear_framework        the nuclear framework that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const NuclearFramework& nuclear_framework);



    // PUBLIC METHODS

    /**
     *  @return the sum of all the charges of the nuclei
     */
    size_t totalNucleicCharge() const;

    /**
     *  @return the nuclei in this nuclear framework as a std::vector
     */
    const std::vector<Nucleus>& nucleiAsVector() const { return this->nuclei; }

    /**
     *  @return the number of nuclei in this nuclear framework
     */
    size_t numberOfNuclei() const { return this->nucleiAsVector().size(); }

    /**
     *  @param index1   the index of the first nucleus
     *  @param index2   the index of the second nucleus
     *
     *  @return the distance between the two nuclei at index1 and index2 in bohr
     */
    double internuclearDistance(const size_t index1, const size_t index2) const;
};


}  // namespace GQCP



#endif  // GQCP_NUCLEARFRAMEWORK_HPP
