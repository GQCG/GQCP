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

#include "Basis/ScalarBasis/GTOBasisSet.hpp"

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param basisset_name                the name of the basisset
 */
GTOBasisSet::GTOBasisSet(const std::string& basisset_name) :
    basisset_name {basisset_name} {}


/**
 *  PUBLIC METHODS
 */

/**
 *  @param nuclear_framework            the nuclear framework containing the nuclei on which the shells should be centered
 * 
 *  @return the shell set by placing the shells corresponding to the basisset information on every nucleus of the nuclear framework
 */
ShellSet<GTOShell> GTOBasisSet::generate(const NuclearFramework& nuclear_framework) const {

    const auto& nuclei = nuclear_framework.nucleiAsVector();
    libint2::BasisSet libint_basis {basisset_name, LibintInterfacer::get().interface(nuclei)};

    return LibintInterfacer::get().interface(libint_basis, nuclei);
}


/**
 *  @param molecule             the molecule containing the nuclei on which the shells should be centered
 * 
 *  @return the shell set by placing the shells corresponding to the basisset information on every nucleus of the molecule
 */
ShellSet<GTOShell> GTOBasisSet::generate(const Molecule& molecule) const {

    return this->generate(molecule.nuclearFramework());
}


}  // namespace GQCP
