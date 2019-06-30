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
#ifndef GQCP_BASENUCLEAROPERATOR_HPP
#define GQCP_BASENUCLEAROPERATOR_HPP


#include "Atom.hpp"

#include <vector>



namespace GQCP {


/**
 *  A base class to represent nuclear first-quantized operators
 */
class BaseNuclearOperator {
protected:
    std::vector<Atom> atoms;  // the atoms that represent the nuclear framework


public:
    // CONSTRUCTORS

    /**
     *  @param atoms                the atoms that represent the nuclear framework
     */
    BaseNuclearOperator(const std::vector<Atom>& atoms);


    // DESTRUCTOR
    virtual ~BaseNuclearOperator() = 0;  // provide a pure virtual destructor to make the class abstract


    // PUBLIC METHODS

    /**
     *  @return the nuclear framework upon which this operator is built
     */
    const std::vector<Atom>& nuclearFramework() const { return this->atoms; }
};


}  // namespace GQCP



#endif  // GQCP_BASENUCLEAROPERATOR_HPP
