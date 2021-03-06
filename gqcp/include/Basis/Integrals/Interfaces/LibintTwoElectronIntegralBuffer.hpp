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


namespace GQCP {


/**
 *  A buffer for storing libint two-electron integrals
 * 
 *  @tparam _N              the number of components the operator has
 */
template <size_t _N>
class LibintTwoElectronIntegralBuffer: public BaseTwoElectronIntegralBuffer<double, _N> {
public:
    using IntegralScalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;   // the number of components the operator has


private:
    using libint2_buffer_t = LibintInterfacer::libint_target_ptr_vec;  // t for type
    libint2_buffer_t libint2_buffer;


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param libint2_buffer       the libint2 buffer that contains the calculated integrals
     *  @param nbf1                 the number of basis functions in the first shell
     *  @param nbf2                 the number of basis functions in the second shell
     *  @param nbf3                 the number of basis functions in the third shell
     *  @param nbf4                 the number of basis functions in the fourth shell
     */
    LibintTwoElectronIntegralBuffer(const libint2_buffer_t& libint2_buffer, const size_t nbf1, const size_t nbf2, const size_t nbf3, const size_t nbf4) :
        libint2_buffer {libint2_buffer},
        BaseTwoElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2, nbf3, nbf4) {}


    /**
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return if all the values of the calculated integrals are zero
     */
    bool areIntegralsAllZero() const override { return (this->libint2_buffer[0] == nullptr); }

    /**
     *  @param i            the index of the component of the operator
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     *  @param f3           the index of the basis function within shell 3
     *  @param f4           the index of the basis function within shell 4
     * 
     *  @return a value from this integral buffer
     */
    IntegralScalar value(const size_t i, const size_t f1, const size_t f2, const size_t f3, const size_t f4) const override {

        return this->libint2_buffer[i][f4 + this->nbf4 * (f3 + this->nbf3 * (f2 + this->nbf2 * (f1)))];  // integrals are packed in row-major form
    }
};


}  // namespace GQCP
