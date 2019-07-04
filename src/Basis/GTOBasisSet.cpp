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
#include "Basis/GTOBasisSet.hpp"

#include "Basis/LibintInterfacer.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param basisset_name                the name of the basisset
 */
GTOBasisSet::GTOBasisSet(const std::string& basisset_name) :
    basisset_name (basisset_name)
{}



/**
 *  PUBLIC METHODS
 */

/**
 *  @param molecule             the molecule containing the nuclei on which the shells should be centered
 * 
 *  @return the shell set by placing the shells corresponding to the basisset information on every nucleus of the molecule
 */
ShellSet<GTOShell> GTOBasisSet::generate(const Molecule& molecule) const {

    // TODO no longer use libint2 to read this

    const auto& nuclei = molecule.nuclearFramework().nucleiAsVector();

    libint2::BasisSet libint_basis (basisset_name, LibintInterfacer::get().interface(nuclei));

    return LibintInterfacer::get().interface(libint_basis, nuclei);
}


}  // namespace GQCP
