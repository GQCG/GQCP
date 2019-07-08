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
#ifndef GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP
#define GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP


#include "Basis/BaseOneElectronIntegralBuffer.hpp"

#include "Basis/LibintInterfacer.hpp"


namespace GQCP {


/**
 *  A buffer for storing libint one-electron integrals
 * 
 *  @tparam _N              the number of components the operator has
 */
template <size_t _N>
class LibintOneElectronIntegralBuffer : public BaseOneElectronIntegralBuffer<double, _N> {
public:
    using IntegralScalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has


private:
    size_t component_offset;  // the number of libint components that should be skipped during access of calculated values (libint2::Operator::emultipole1 has 4 libint2 components, but in reality there should only be 1)
    double scaling_factor;  // a factor that is multiplied to all of the calculated integrals

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
     *  @param 
     */
    LibintOneElectronIntegralBuffer(const libint2_buffer_t& libint2_buffer, const size_t nbf1, const size_t nbf2, const size_t component_offset=0, const double scaling_factor=1.0) :
        libint2_buffer (libint2_buffer),
        component_offset (component_offset),
        scaling_factor (scaling_factor),
        BaseOneElectronIntegralBuffer<IntegralScalar, N>(nbf1, nbf2)
    {}



    /**
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return the matrix representation of the integrals that are in this buffer
     */
    std::array<Matrix<IntegralScalar>, N> integrals() const override {

        // Initialize N zero partial components
        std::array<Matrix<IntegralScalar>, N> partial_components;  // N partial components of the total matrix representation of the operator
        for (auto& partial_component : partial_components) {
            partial_component = Matrix<IntegralScalar>::Zero(this->nbf1, this->nbf2);
        }


        // Place the calculated integrals inside the partial components
        for (size_t f1 = 0; f1 != this->nbf1; f1++) {  // f1: index of basis function within shell 1
            for (size_t f2 = 0; f2 != this->nbf2; f2++) {  // f2: index of basis function within shell 2

                for (size_t i = 0; i < N; i++) {
                    partial_components[i](f1, f2) = this->extract(i, f1, f2);
                }

            }
        }  // data access loops

        return partial_components;
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param i            the operator component number
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     * 
     *  @return the integral corresponding to the given basis functions from the buffer
     */
    const IntegralScalar& extract(const size_t i, const size_t f1, const size_t f2) const {
        return this->scaling_factor * this->libint2_buffer[i + this->component_offset][f2 + f1 * this->nbf2]; // integrals are packed in row-major form
    }
};


}  // namespace GQCP



#endif  // GQCP_LIBINTONEELECTRONINTEGRALBUFFER_HPP
