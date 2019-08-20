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
#ifndef GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP
#define GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP


#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"

#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"


namespace GQCP {


/**
 *  An one-electron integral engine that uses libcint as its backend
 * 
 *  @tparam _ShellType                  the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  _ShellType is a template parameter because that enables compile-time checking of correct arguments
 */
template <typename _ShellType, size_t _N, typename _IntegralScalar>
class LibcintOneElectronIntegralEngine : public BaseOneElectronIntegralEngine<_ShellType, _N, _IntegralScalar> {
public:
    using ShellType = _ShellType;  // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


private:
    Libcint1eFunction libcint_function;  // the libcint one-electron integral function

    // Parameters to pass to the buffer
    double scaling_factor = 1.0;  // a factor that is multiplied to all of the calculated integrals


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param op               the overlap operator
     */
    LibcintOneElectronIntegralEngine(const OverlapOperator& op) :
        libcint_function (LibcintInterfacer().oneElectronFunction(op))
    {}

    /**
     *  @param op               the kinetic operator
     */
    LibcintOneElectronIntegralEngine(const KineticOperator& op) :
        libcint_function (LibcintInterfacer().oneElectronFunction(op))
    {}

    /**
     *  @param op               the nuclear attraction operator
     */
    LibcintOneElectronIntegralEngine(const NuclearAttractionOperator& op) :
        libcint_function (LibcintInterfacer().oneElectronFunction(op))
    {}

    /**
     *  @param op               the electronic electric dipole operator
     */
    LibcintOneElectronIntegralEngine(const ElectronicDipoleOperator& op) :
        libcint_function (LibcintInterfacer().oneElectronFunction(op)),
        scaling_factor (-1.0)  // apply the minus sign which comes from the charge of the electrons -e
    {}



    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @param shell1           the first shell
     *  @param shell2           the second shell
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const GTOShell& shell1, const GTOShell& shell2) override {

        // Interface the given GTOShells to libcint-compatible data
        const LibcintInterfacer libcint_interfacer;
        auto libcint_raw_container = libcint_interfacer.convert({shell1, shell2});
        int shell_indices[2] = {0, 1};  // we're interfacing for every shell pair, so the shell indices are 0 and 1


        // Pre-allocate a raw buffer, because libcint functions expect a data pointer
        const size_t nbf1 = shell1.numberOfBasisFunctions();
        const size_t nbf2 = shell2.numberOfBasisFunctions();
        double libcint_buffer[N * nbf1 * nbf2];


        // Let libcint compute the integrals and return the corresponding buffer
        this->libcint_function(libcint_buffer, shell_indices, libcint_raw_container.atmData(), libcint_raw_container.numberOfAtoms(), libcint_raw_container.basData(), libcint_raw_container.numberOfBasisFunctions(), libcint_raw_container.envData());
        std::vector<double> buffer_converted (libcint_buffer, libcint_buffer + N*nbf1*nbf2);  // std::vector constructor from .begin() and .end()

        return std::make_shared<LibcintOneElectronIntegralBuffer<IntegralScalar, N>>(buffer_converted, nbf1, nbf2, this->scaling_factor);
    }
};



}  // namespace GQCP



#endif  // GQCP_LIBCINTONEELECTRONINTEGRALENGINE_HPP
