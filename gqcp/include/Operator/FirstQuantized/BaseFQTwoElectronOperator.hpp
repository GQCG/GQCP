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


#include <cstddef>


namespace GQCP {


/**
 *  A base class used to represent two-electron operators.
 * 
 *  @tparam _Scalar             the scalar representation of the operator
 *  @tparam _Components         the number of components the operator has
 */
template <typename _Scalar, size_t _Components>
class BaseFQTwoElectronOperator {
public:
    using Scalar = _Scalar;                          // the scalar representation of the operator
    static constexpr auto Components = _Components;  // the number of components the operator has


public:
    // DESTRUCTORS

    virtual ~BaseFQTwoElectronOperator() {}
};


/*
 *  CONVENIENCE ALIASES
 */
template <typename Scalar>
using BaseScalarFQTwoElectronOperator = BaseFQTwoElectronOperator<Scalar, 1>;

template <typename Scalar>
using BaseVectorFQTwoElectronOperator = BaseFQTwoElectronOperator<Scalar, 3>;


}  // namespace GQCP
