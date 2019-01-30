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
#include "AOBasis.hpp"
#include "LibintCommunicator.hpp"

namespace GQCP {


/**
 *  Constructor that creates an AOBasis from a basisset name
 *
 *  @param molecule         the molecule to which the AO basis corresponds
 *  @param basis_set        the name of the basisset
 */
AOBasis::AOBasis(const Molecule& molecule, const std::string& basis_set) :
    atoms (molecule.get_atoms()),
    basis_functions (libint2::BasisSet(std::move(basis_set), LibintCommunicator::get().interface(this->atoms))),  // construct a libint2::BasisSet
    number_of_basis_functions (static_cast<size_t>(this->basis_functions.nbf()))
{}


}  // namespace GQCP
