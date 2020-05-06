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

#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param nuclear_framework            the nuclear framework underlying a nuclear operator
 */
BaseNuclearOperator::BaseNuclearOperator(const NuclearFramework& nuclear_framework) :
    nuclear_framework {nuclear_framework} {}


/**
 *  @param nuclei                       the nuclei that are considered to represent the nuclear framework underlying a nuclear operator
 */
BaseNuclearOperator::BaseNuclearOperator(const std::vector<Nucleus>& nuclei) :
    BaseNuclearOperator(NuclearFramework(nuclei)) {}


/*
 *  DESTRUCTOR
 */

BaseNuclearOperator::~BaseNuclearOperator() {};  // pure virtual destructor should have an empty implementation


}  // namespace GQCP
