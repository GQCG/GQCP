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


#include "Molecule/NuclearFramework.hpp"


namespace GQCP {


/**
 *  A base class to represent nuclear first-quantized operators
 */
class BaseNuclearOperator {
protected:
    NuclearFramework nuclear_framework;  // the nuclear framework underlying a nuclear operator


public:
    // CONSTRUCTORS

    /**
     *  @param nuclear_framework            the nuclear framework underlying a nuclear operator
     */
    BaseNuclearOperator(const NuclearFramework& nuclear_framework);

    /**
     *  @param nuclei                       the nuclei that are considered to represent the nuclear framework underlying a nuclear operator
     */
    BaseNuclearOperator(const std::vector<Nucleus>& nuclei);


    // DESTRUCTOR
    virtual ~BaseNuclearOperator() = 0;  // provide a pure virtual destructor to make the class abstract


    // PUBLIC METHODS

    /**
     *  @return the nuclear framework upon which this operator is built
     */
    const NuclearFramework& nuclearFramework() const { return this->nuclear_framework; }
};


}  // namespace GQCP
