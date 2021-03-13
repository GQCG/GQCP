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


#include "Basis/Integrals/BaseTwoElectronIntegralBuffer.hpp"

#include <array>


namespace GQCP {


/**
 *  A buffer for storing two-electron integrals.
 * 
 *  An integral buffer facilitates placing integrals (see emplace()) that were calculated over shells into the correct matrix representation of an operator in an orbital basis.
 * 
 *  @tparam _IntegralScalar         The scalar representation of an integral.
 *  @tparam _N                      The number of components the corresponding operator has.
 */
template <typename _IntegralScalar, size_t _N>
class TwoElectronIntegralBuffer:
    public BaseTwoElectronIntegralBuffer<_IntegralScalar, _N> {
public:
    // The scalar representation of an integral.
    using IntegralScalar = _IntegralScalar;

    // The number of components the corresponding operator has.
    static constexpr auto N = _N;


protected:
    // The calculated integrals.
    std::array<std::vector<IntegralScalar>, N> buffer;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param nbf1             The number of basis functions in the first shell.
     *  @param nbf2             The number of basis functions in the second shell.
     *  @param nbf3             The number of basis functions in the third shell.
     *  @param nbf4             The number of basis functions in the fourth shell.
     *  @param buffer           The calculated integrals.
     */
    TwoElectronIntegralBuffer(const size_t nbf1, const size_t nbf2, const size_t nbf3, const size_t nbf4, const std::array<std::vector<IntegralScalar>, N>& buffer) :
        buffer {buffer},
        BaseTwoElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2, nbf3, nbf4) {}


    /*
     *  MARK: Emplacement
     */

    /**
     *  @return If all the values of the calculated integrals are zero.
     */
    bool areIntegralsAllZero() const { return false; }

    /**
     *  @param i            The index of the component of the operator.
     *  @param f1           The index of the basis function within shell 1.
     *  @param f2           The index of the basis function within shell 2.
     *  @param f3           The index of the basis function within shell 3.
     *  @param f4           The index of the basis function within shell 4.
     * 
     *  @return A value from this integral buffer.
     */
    IntegralScalar value(const size_t i, const size_t f1, const size_t f2, const size_t f3, const size_t f4) const override {
        // We'll access the component first, then a row-major ordering of the calculated integrals.
        return this->buffer[i][f4 + this->nbf4 * (f3 + this->nbf3 * (f2 + this->nbf2 * f1))];
    }
};


}  // namespace GQCP
