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


#include <cstddef>


namespace GQCP {


/**
 *  A base class used to represent one-electron operators
 * 
 *  @tparam _Scalar         the scalar representation of the operator
 *  @tparam _Components     the number of components the operator has
 */
template <typename _Scalar, size_t _Components>
class BaseFQOneElectronOperator {
public:
    using Scalar = _Scalar;  // the scalar representation of the operator
    static constexpr auto Components = _Components;  // the number of components the operator has


public:
    // DESTRUCTOR

    virtual ~BaseFQOneElectronOperator() {}
};



/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using BaseScalarFQOneElectronOperator = BaseFQOneElectronOperator<Scalar, 1>;

template <typename Scalar>
using BaseVectorFQOneElectronOperator = BaseFQOneElectronOperator<Scalar, 3>;


}  // namespace GQCP
