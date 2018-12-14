// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_AOBASIS_HPP
#define GQCP_AOBASIS_HPP


#include "Atom.hpp"
#include "Molecule.hpp"

#include <Eigen/Dense>
#include <libint2.hpp>




namespace GQCP {


/**
 *  A class that represents an atomic orbital basis
 */
class AOBasis {
private:
    std::vector<Atom> atoms;
    libint2::BasisSet basis_functions;
    size_t number_of_basis_functions;

public:
    // CONSTRUCTOR
    /**
     *  Constructor that creates an AOBasis from a basisset name
     *
     *  @param molecule         the molecule to which the AO basis corresponds
     *  @param basis_set        the name of the basisset
     */
    AOBasis(const Molecule& molecule, const std::string& basis_set);


    // GETTERS
    const std::vector<Atom>& get_atoms() const { return this->atoms; }
    const libint2::BasisSet& get_basis_functions() const { return this->basis_functions; }
    size_t get_number_of_basis_functions() const { return this->number_of_basis_functions; }
};


}  // namespace GQCP


#endif  // GQCP_AOBASIS_HPP
