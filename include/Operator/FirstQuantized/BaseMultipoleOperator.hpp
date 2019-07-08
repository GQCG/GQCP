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
#ifndef GQCP_BASEMULTIPOLEOPERATOR_HPP
#define GQCP_BASEMULTIPOLEOPERATOR_HPP


#include "Mathematical/Matrix.hpp"



namespace GQCP {


/**
 *  A base class to represent first-quantized multipole operators
 */
class BaseMultipoleOperator {
protected:
    Vector<double, 3> o;  // the origin of the multipole


public:
    // CONSTRUCTORS

    /**
     *  @param o        the origin of the multipole
     */
    BaseMultipoleOperator(const Vector<double, 3>& o=Vector<double, 3>::Zero());


    // DESTRUCTOR
    virtual ~BaseMultipoleOperator() = 0;


    // PUBLIC METHODS

    /**
     *  @return the origin of the multipole operator
     */
    const Vector<double, 3>& origin() const { return this->o; }
};


}  // namespace GQCP



#endif  // GQCP_BASEMULTIPOLEOPERATOR_HPP
