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

#include <array>


namespace GQCP {


/**
 *  A base class to implement buffers for storing one-electron integrals
 * 
 *  An integral buffer facilitates placing integrals (see emplace()) that were calculated over shells into the correct matrix representation of an operator in an orbital basis
 * 
 *  @tparam _IntegralScalar         the scalar representation of an integral
 *  @tparam _N                      the number of components the corresponding operator has
 */
template <typename _IntegralScalar, size_t _N>
class OneElectronIntegralBuffer:
    public BaseOneElectronIntegralBuffer<_IntegralScalar, _N> {
public:
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;            // the number of components the operator has


protected:
    std::array<std::vector<IntegralScalar>, N> buffer;  // the calculated integrals

public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param nbf1             the number of basis functions in the first shell
     *  @param nbf2             the number of basis functions in the second shell
     *  @param buffer           the calculated integrals
     */
    OneElectronIntegralBuffer(const size_t nbf1, const size_t nbf2, const std::array<std::vector<IntegralScalar>, N>& buffer) :
        buffer {buffer},
        BaseOneElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2) {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return if all the values of the calculated integrals are zero
     */
    bool areIntegralsAllZero() const { return false; }

    /**
     *  @param i            the index of the component of the operator
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     * 
     *  @return a value from this integral buffer
     */
    virtual IntegralScalar value(const size_t i, const size_t f1, const size_t f2) const {
        return this->buffer[i][f2 + f1 * this->nbf2];  // accessing the component first, then row-major ordering of the calculated integrals
    }
};


}  // namespace GQCP
