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
#ifndef GQCP_BASEHAMILTONIANPARAMETERS_HPP
#define GQCP_BASEHAMILTONIANPARAMETERS_HPP


#include <memory>

#include "AOBasis.hpp"


namespace GQCP {


class BaseHamiltonianParameters {
protected:
    std::shared_ptr<GQCP::AOBasis> ao_basis;  // the initial atomic orbitals

public:
    // CONSTRUCTOR
    /**
     *  Constructor based on a given @param ao_basis
     */
    explicit BaseHamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseHamiltonianParameters() = 0;
};


}  // namespace GQCP


#endif  // GQCP_BASEHAMILTONIANPARAMETERS_HPP
