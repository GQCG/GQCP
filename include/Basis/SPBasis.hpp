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


namespace GQCP {


/**
 *  A class that represents a single-particle basis
 * 
 *  @tparam _SBasis                     the scalar basis that underlies the single-particle basis
 *  @tparam _TransformationScalar       the scalar type of the transformation matrix that connects the scalar basis with the current single-particle 'orbitals'
 */
template <typename _SBasis, typename _TransformationScalar>
class SPBasis {
public:
    using SBasis = _SBasis;  // can't use ScalarBasis because that's the name of the class
    using TransformationScalar = _TransformationScalar;


private:



public:
    /*
     *  PUBLIC METHODS
     */

    
};


}  // namespace GQCP
