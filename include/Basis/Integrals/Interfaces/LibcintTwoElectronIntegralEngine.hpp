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
#ifndef GQCP_LIBCINTTWOELECTRONINTEGRALENGINE_HPP
#define GQCP_LIBCINTTWOELECTRONINTEGRALENGINE_HPP


#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"

#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"


namespace GQCP {


/**
 *  An two-electron integral engine that uses libcint as its backend
 * 
 *  @tparam _ShellType                  the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  _ShellType is a template parameter because that enables compile-time checking of correct arguments
 */
template <typename _ShellType, size_t _N, typename _IntegralScalar>
class LibcintTwoElectronIntegralEngine : public BaseTwoElectronIntegralEngine<_ShellType, _N, _IntegralScalar> {
public:
    using ShellType = _ShellType;  // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


private:
    Libcint2eFunction libcint_function;  // the libcint two-electron integral function


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param op           the Coulomb repulsion operator
     */
    LibcintTwoElectronIntegralEngine(const CoulombRepulsionOperator& op) :
        libcint_function (LibcintInterfacer().twoElectronFunction(op))
    {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @param shell1          the first shell
     *  @param shell2          the second shell
     *  @param shell3          the third shell
     *  @param shell4          the fourth shell
     * 
     *  This method is not marked const to allow the Engine's internals to be changed
     */
    std::shared_ptr<BaseTwoElectronIntegralBuffer<IntegralScalar, N>> calculate(const ShellType& shell1, const ShellType& shell2, const ShellType& shell3, const ShellType& shell4) override {

        // Interface the given GTOShells to libcint-compatible data
        const LibcintInterfacer libcint_interfacer;
        auto libcint_raw_container = libcint_interfacer.convert({shell1, shell2, shell3, shell4});
        int shell_indices[4] = {0, 1, 2, 3};  // we're interfacing for every shell quartet, so the shell indices are 0, 1, 2, and 3


        // Pre-allocate a raw buffer, because libcint functions expect a data pointer
        const size_t nbf1 = shell1.numberOfBasisFunctions();
        const size_t nbf2 = shell2.numberOfBasisFunctions();
        const size_t nbf3 = shell3.numberOfBasisFunctions();
        const size_t nbf4 = shell4.numberOfBasisFunctions();
        double libcint_buffer[N * nbf1 * nbf2 * nbf3 * nbf4];


        // Let libcint compute the integrals and return the corresponding buffer
        this->libcint_function(libcint_buffer, shell_indices, libcint_raw_container.atmData(), libcint_raw_container.numberOfAtoms(), libcint_raw_container.basData(), libcint_raw_container.numberOfBasisFunctions(), libcint_raw_container.envData(), nullptr);
        std::vector<double> buffer_converted (libcint_buffer, libcint_buffer + N*nbf1*nbf2*nbf3*nbf4);  // std::vector constructor from .begin() and .end()

        return std::make_shared<LibcintTwoElectronIntegralBuffer<IntegralScalar, N>>(buffer_converted, nbf1, nbf2, nbf3, nbf4);
    }
};


}  // namespace GQCP



#endif  // GQCP_LIBCINTTWOELECTRONINTEGRALENGINE_HPP
