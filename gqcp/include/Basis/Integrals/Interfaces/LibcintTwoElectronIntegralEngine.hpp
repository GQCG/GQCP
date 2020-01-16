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
#pragma once


#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"

#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralBuffer.hpp"
#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/**
 *  An two-electron integral engine that uses libcint as its backend
 * 
 *  @tparam _Shell                      the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  @note _Shell is a template parameter because that enables compile-time checking of correct arguments.
 *  See also the notes in LibcintOneElectronIntegralEngine.
 *  The libcint optimizer struct should also be kept in the engine during the shell-quartet loop, because it should only be initialized once, having access to all the data inside the libcint RawContainer.
 */
template <typename _Shell, size_t _N, typename _IntegralScalar>
class LibcintTwoElectronIntegralEngine : public BaseTwoElectronIntegralEngine<_Shell, _N, _IntegralScalar> {
public:
    using Shell = _Shell;  // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


private:
    Libcint2eFunction libcint_function;  // the libcint two-electron integral function
    Libcint2eOptimizerFunction libcint_optimizer_function;  // the libcint two-electron optimizer integral function

    // Data that has to be kept as a member (see the class note)
    libcint::RawContainer libcint_raw_container;  // the raw libcint data
    ShellSet<Shell> shell_set;  // the corresponding shell set


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     */
    LibcintTwoElectronIntegralEngine(const CoulombRepulsionOperator& op, const ShellSet<Shell>& shell_set) :
        libcint_function (LibcintInterfacer().twoElectronFunction(op)),
        libcint_optimizer_function (LibcintInterfacer().twoElectronOptimizerFunction(op)),
        libcint_raw_container (LibcintInterfacer().convert(shell_set)),
        shell_set (shell_set)
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
    std::shared_ptr<BaseTwoElectronIntegralBuffer<IntegralScalar, N>> calculate(const Shell& shell1, const Shell& shell2, const Shell& shell3, const Shell& shell4) override {

        // Find to which indices in the RawContainer the given shells correspond
        int shell_indices[4];
        shell_indices[0] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell1));
        shell_indices[1] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell2));
        shell_indices[2] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell3));
        shell_indices[3] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell4));


        // Pre-allocate a raw buffer, because libcint functions expect a data pointer
        const size_t nbf1 = shell1.numberOfBasisFunctions();
        const size_t nbf2 = shell2.numberOfBasisFunctions();
        const size_t nbf3 = shell3.numberOfBasisFunctions();
        const size_t nbf4 = shell4.numberOfBasisFunctions();
        double libcint_buffer[N * nbf1 * nbf2 * nbf3 * nbf4];


        // Let libcint compute the integrals and return the corresponding buffer
        const auto result = this->libcint_function(libcint_buffer, shell_indices, this->libcint_raw_container.atmData(), this->libcint_raw_container.numberOfAtoms(), this->libcint_raw_container.basData(), this->libcint_raw_container.numberOfBasisFunctions(), this->libcint_raw_container.envData(), nullptr);  // no optimizer struct


        std::vector<double> buffer_converted (libcint_buffer, libcint_buffer + N*nbf1*nbf2*nbf3*nbf4);  // std::vector constructor from .begin() and .end()
        return std::make_shared<LibcintTwoElectronIntegralBuffer<IntegralScalar, N>>(buffer_converted, nbf1, nbf2, nbf3, nbf4, result);
    }
};


}  // namespace GQCP
