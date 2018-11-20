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
#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param ao_basis     the initial AO basis
 *  @param scalar       the scalar interaction term
 */
BaseHamiltonianParameters::BaseHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis, double scalar) :
    ao_basis (std::move(ao_basis)),
    scalar (scalar)
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseHamiltonianParameters::~BaseHamiltonianParameters() {}



}  // namespace GQCP
