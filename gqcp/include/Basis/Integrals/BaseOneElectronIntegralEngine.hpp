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


#include "Basis/Integrals/BaseOneElectronIntegralBuffer.hpp"

#include <memory>


namespace GQCP {


/**
 *  A base class that helps to implement one-electron integral engines. Integral engines are used calculate integrals of operators over shells. See also the calculate() call.
 * 
 *  @tparam _Shell                  The type of shell the integral engine is able to handle. This enables compile-time checking of correct arguments.
 *  @tparam _N                      The number of components the operator has.
 *  @tparam _IntegralScalar         
 *  */
template <typename _Shell, size_t _N, typename _IntegralScalar>
class BaseOneElectronIntegralEngine {
public:
    // The type of shell the integral engine is able to handle. This enables compile-time checking of correct arguments.
    using Shell = _Shell;

    // The number of components the operator has.
    static constexpr auto N = _N;  // the number of components the operator has

    // The scalar representation of an integral.
    using IntegralScalar = _IntegralScalar;


public:
    /**
     *  A default destructor.
     */
    virtual ~BaseOneElectronIntegralEngine() = default;


    /*
     *  MARK: Integral calculations
     */

    /**
     *  Calculate all the integrals over the given shells.
     * 
     *  @param shell1           The first shell.
     *  @param shell2           The second shell.
     * 
     *  @return A buffer containing the calculated integrals.
     * 
     *  @note This method is not marked const to allow the Engine's internals to be changed.
     */
    virtual std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const Shell& shell1, const Shell& shell2) = 0;
};


}  // namespace GQCP
