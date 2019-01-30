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
#ifndef GQCP_BASERDM_HPP
#define GQCP_BASERDM_HPP


#include <cstdlib>


namespace GQCP {


/**
 *  A base class for the representation of reduced density matrices
 */
class BaseRDM {
protected:
    size_t dim;  // dimension of the matrix representation of the operator


public:
    // CONSTRUCTORS
    /**
     *  @param dimension    the dimension of the matrix representation of the operator (i.e. the number of orbitals)
     */
    explicit BaseRDM(size_t dimension);


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseRDM() = 0;


    // GETTERS
    size_t get_dim() const { return this->dim; }
};


}  // namespace GQCP


#endif  // GQCP_BASERDM_HPP
