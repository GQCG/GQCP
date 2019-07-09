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
#ifndef GQCP_BASETWOELECTRONINTEGRALENGINE_HPP
#define GQCP_BASETWOELECTRONINTEGRALENGINE_HPP


#include "Basis/BaseTwoElectronIntegralBuffer.hpp"


namespace GQCP {


/**
 *  A base class to implement two-electron integral engines. Integral engines are used calculate integrals of operators over shells, see also the calculate() call
 * 
 *  @tparam _ShellType                  the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  _ShellType is a template parameter because that enables compile-time checking of correct arguments
 */
template <typename _ShellType, size_t _N, typename _IntegralScalar>
class BaseTwoElectronIntegralEngine {
public:
    using ShellType = _ShellType;  // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


public:
    // PUBLIC METHODS

    /**
     *  @param sh1          the first shell
     *  @param sh2          the second shell
     *  @param sh3          the third shell
     *  @param sh4          the fourth shell
     * 
     *  This method is not marked const to allow the Engine's internals to be changed
     */
    virtual std::shared_ptr<BaseTwoElectronIntegralBuffer<IntegralScalar, N>> calculate(const ShellType& sh1, const ShellType& sh2, const ShellType& sh3, const ShellType& sh4) = 0;
};


}  // namespace GQCP



#endif  // GQCP_BASETWOELECTRONINTEGRALENGINE_HPP
