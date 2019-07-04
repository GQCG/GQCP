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
#ifndef GQCP_GTOBASISSETFACTORY_HPP
#define GQCP_GTOBASISSETFACTORY_HPP


#include "Basis/GTOShell.hpp"
#include "Basis/ShellSet.hpp"
#include "Molecule/Molecule.hpp"

#include <string>


namespace GQCP {



/**
 *  A class that represents a Gaussian basis set: it serves as a mold to construct GTOShells
 */
class GTOBasisSet {
private:
    std::string basisset_name;  // the name of the basisset


public:
    // CONSTRUCTORS

    /**
     *  @param basisset_name                the name of the basisset
     */
    GTOBasisSet(const std::string& basisset_name);


    // PUBLIC METHODS

    /**
     *  @return the name of the basisset
     */
    const std::string& name() const { return this->basisset_name; }

    /**
     *  @param molecule             the molecule containing the nuclei on which the shells should be centered
     * 
     *  @return the shell set by placing the shells corresponding to the basisset information on every nucleus of the molecule
     */
    ShellSet<GTOShell> generate(const Molecule& molecule) const;
};


}  // namespace GQCP



#endif  // GQCP_GTOBASISSETFACTORY_HPP
