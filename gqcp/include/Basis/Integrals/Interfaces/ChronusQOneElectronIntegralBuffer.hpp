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
#include "Utilities/aliases.hpp"

namespace GQCP {


/**
 *  A buffer for storing one-electron integrals produced by a ChronusQ engine.
 * 
 *  @tparam _N              The number of components the operator has.
 */
template <size_t _N>
class ChronusQOneElectronIntegralBuffer:
    public BaseOneElectronIntegralBuffer<complex, _N> {
public:
    // The type used to represent an integral calculated by ChronusQ.
    using IntegralScalar = complex;

    // The number of components the operator has.
    static constexpr auto N = _N;


private:
    // The ChronusQ-calculated integrals over a pair of shells.
    using chronusq_buffer_t = std::vector<std::vector<complex>>;  // TODO: Move to ChronusQInterfacer
    chronusq_buffer_t chronusq_buffer;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Initialize a `ChronusQOneElectronIntegralBuffer` from a set of calculated integrals and numbers of basis functions.
     * 
     *  @param chronusq_buffer          The ChronusQ-calculated integrals over a pair of shells.
     *  @param nbf1                     The number of basis functions in the first shell.
     *  @param nbf2                     The number of basis functions in the second shell.
     */
    ChronusQOneElectronIntegralBuffer(const chronusq_buffer_t& chronusq_buffer, const size_t nbf1, const size_t nbf2) :
        chronusq_buffer {chronusq_buffer},
        BaseOneElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2) {}


    /**
     *  MARK: Access
     */

    /**
     *  @return If all the values of the calculated integrals are zero.
     */
    bool areIntegralsAllZero() const override { return false; /* Don't let ChronusQ assume any integrals are zero. */ }

    /**
     *  @param i            The index of the component.
     *  @param f1           The index of the basis function within shell 1.
     *  @param f2           The index of the basis function within shell 2.
     * 
     *  @return The integral over the given basis functions.
     */
    IntegralScalar value(const size_t i, const size_t f1, const size_t f2) const override {
        return this->chronusq_buffer[i][f2 + f1 * this->nbf2];  // Like Libint2, ChronusQ integrals are packed in row-major form.
    }
};


}  // namespace GQCP
