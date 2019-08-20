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
#ifndef GQCP_LIBCINTONEELECTRONINTEGRALBUFFER_HPP
#define GQCP_LIBCINTONEELECTRONINTEGRALBUFFER_HPP


#include "Basis/Integrals/BaseOneElectronIntegralBuffer.hpp"



namespace GQCP {


/**
 *  A buffer for storing Libcint one-electron integrals
 * 
 *  @tparam _IntegralScalar         the scalar representation of an integral
 *  @tparam _N                      the number of components the operator has
 * 
 *  @note Since a Libcint function accepts a raw pointer/buffer as an argument, this class should take ownership of the buffer data since the raw buffer is created on the stack. This is done by placing the raw buffer data inside a std::vector.
 */
template <typename _IntegralScalar, size_t _N>
class LibcintOneElectronIntegralBuffer : public BaseOneElectronIntegralBuffer<_IntegralScalar, _N> {
public:
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


private:
    double scaling_factor;  // a factor that multiplies every calculated value (for example in the dipole integrals, an extra factor -1 should be included)
    std::vector<IntegralScalar> buffer;  // the libcint integral data converted to a C++ vector


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param buffer               the libcint integral data, already put inside a vector
     *  @param nbf1                 the number of basis functions in the first shell
     *  @param nbf2                 the number of basis functions in the second shell
     *  @param scaling_factor       a factor that multiplies every calculated value (for example in the dipole integrals, an extra factor -1 should be included)
     */
    LibcintOneElectronIntegralBuffer(const std::vector<IntegralScalar>& buffer, const size_t nbf1, const size_t nbf2, const double scaling_factor=1.0) :
        scaling_factor (scaling_factor),
        buffer (buffer),
        BaseOneElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2)
    {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */
    
    /**
     *  @param i            the index of the component of the operator
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     * 
     *  @return a value from this integral buffer
     */
    virtual IntegralScalar value(const size_t i, const size_t f1, const size_t f2) const {
        return this->scaling_factor * this->buffer[f1 + this->nbf1 * (f2 + this->nbf2 * i)];
    }
};



}  // namespace GQCP



#endif  // GQCP_LIBCINTONEELECTRONINTEGRALBUFFER_HPP
