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
#ifndef GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP
#define GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP


#include "Basis/BaseOneElectronIntegralEngine.hpp"

#include "Basis/LibcintOneElectronIntegralBuffer.hpp"


namespace GQCP {


/**
 *  An one-electron integral engine that uses libcint as its backend
 * 
 *  @tparam _ShellType                  the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  _ShellType is a template parameter because that enables compile-time checking of correct arguments
 */
template <typename _ShellType, size_t _N, typename _IntegralScalar>
class LibcintOneElectronIntegralEngine : public BaseOneElectronIntegralEngine<_ShellType, _N, _IntegralScalar> {
public:
    using ShellType = _ShellType;  // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has
};



}  // namespace GQCP



#endif  // GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP
