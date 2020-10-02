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

#pragma once


#include "Molecule/NuclearFramework.hpp"


namespace GQCP {


/**
 *  A base class/interface used to represent nuclear first-quantized operators.
 */
class BaseNuclearOperator {
private:
    // The nuclear framework underlying a nuclear operator.
    NuclearFramework nuclear_framework;

public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Construct a `BaseNuclearOperator` from a nuclear framework.
     * 
     *  @param nuclear_framework            The nuclear framework underlying a nuclear operator.
     */
    BaseNuclearOperator(const NuclearFramework& nuclear_framework);

    /**
     *  Construct a `BaseNuclearOperator` from a vector of nuclei.
     * 
     *  @param nuclei                       The nuclei that are considered to represent the nuclear framework underlying a nuclear operator.
     */
    BaseNuclearOperator(const std::vector<Nucleus>& nuclei);


    /*
     *  MARK: Destructor
     */

    // Make the destructor pure virtual in order to make this class abstract.
    virtual ~BaseNuclearOperator() = 0;


    /*
     *  MARK: Nuclear framework
     */

    /**
     *  @return The nuclear framework underlying a nuclear operator.
     */
    const NuclearFramework& nuclearFramework() const { return this->nuclear_framework; }
};


}  // namespace GQCP
