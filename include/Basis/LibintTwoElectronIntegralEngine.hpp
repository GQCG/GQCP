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
// 
#ifndef GQCP_LIBINTTWOELECTRONINTEGRALENGINE_HPP
#define GQCP_LIBINTTWOELECTRONINTEGRALENGINE_HPP


#include "Basis/BaseTwoElectronIntegralEngine.hpp"


namespace GQCP {


/**
 *  A two-electron integral engine that uses libint as its backend
 * 
 *  @tparam _ShellType          the type of shell the integral engine is able to handle
 *  @tparam _N                  the number of components the operator has
 *  @tparam _Scalar             the scalar representation of an integral
 * 
 *  _ShellType is a template parameter because that enables compile-time checking of correct arguments
 */
template <typename _ShellType, size_t _N>
class LibintTwoElectronIntegralEngine : public BaseTwoElectronIntegralEngine<_ShellType, _N, double> {
public:
    using ShellType = _ShellType;  // the type of shell the integral engine is able to handle
    using Scalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has
};


}  // namespace GQCP



#endif  // GQCP_LIBINTTWOELECTRONINTEGRALENGINE_HPP
