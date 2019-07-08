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
#ifndef GQCP_LIBINTTWOELECTRONINTEGRALBUFFER_HPP
#define GQCP_LIBINTTWOELECTRONINTEGRALBUFFER_HPP



#include "Basis/BaseOneElectronIntegralBuffer.hpp"



namespace GQCP {


/**
 *  A buffer for storing libint two-electron integrals
 * 
 *  @tparam _N              the number of components the operator has
 */
template <size_t _N>
class LibintTwoElectronIntegralBuffer : public BaseTwoElectronIntegralBuffer<double, _N> {
public:
    using Scalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has


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
        libint2_buffer (libint2_buffer),
        BaseTwoElectronIntegralBuffer<Scalar, N>(nbf1, nbf2, nbf3, nbf4)
    {}



    /**
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return the matrix representation of the integrals that are in this buffer
     */
    std::array<Tensor<Scalar, 4>, N> integrals() const override {

        // Initialize N zero partial components
        std::array<Tensor<Scalar, 4>, N> partial_components;  // N partial components of the total matrix representation of the operator
        for (auto& partial_component : partial_components) {
            partial_component = Tensor<Scalar, 4>(this->nbf1, this->nbf2, this->nbf3, this->nbf4);
        }


        // Place the calculated integrals inside the partial components
        for (size_t f1 = 0; f1 != this->nbf1; f1++) {  // f1: index of basis function within shell 1
            for (size_t f2 = 0; f2 != this->nbf2; f2++) {  // f2: index of basis function within shell 2
                for (size_t f3 = 0; f3 != this->nbf3; f3++) {  // f3: index of basis function within shell 3
                    for (size_t f4 = 0; f4 != this->nbf4; f4++) {  // f4: index of basis function within shell 4

                        for (size_t i = 0; i < N; i++) {
                            const Scalar calculated_integral = libint2_buffer[i][f4 + this->nbf4 * (f3 + this->nbf3 * (f2 + this->nbf2 * (f1)))];  // integrals are packed in row-major form
                            partial_components[i](f1, f2, f3, f4) = calculated_integral;  // in chemist's notation
                        }

                    }
                }
            }
        }  // data access loops

        return partial_components;
    }
};



}  // namespace GQCP



#endif  // GQCP_LIBINTTWOELECTRONINTEGRALBUFFER_HPP
