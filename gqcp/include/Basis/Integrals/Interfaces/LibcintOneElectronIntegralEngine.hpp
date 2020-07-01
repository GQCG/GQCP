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


#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralBuffer.hpp"
#include "Utilities/literals.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  An one-electron integral engine that uses libcint as its backend
 * 
 *  @tparam _Shell                      the type of shell the integral engine is able to handle
 *  @tparam _N                          the number of components the operator has
 *  @tparam _IntegralScalar             the scalar representation of an integral
 * 
 *  @note _Shell is a template parameter because that enables compile-time checking of correct arguments.
 * 
 *  @note Libcint1eFunction and Libcint2eFunction need the full libcint data (that is stored in a RawContainer), because the electron electrical dipole function needs to access the 'common origin' and the nuclear attraction function needs to know all the nuclei that should be taken into account. We therefore keep the full RawContainer as a member.
 * 
 *  @note In the overridden calculate(const GTOShell&, const GTOShell&) method, the shell indices that correspond to the RawContainer should also be known, which is why we also keep the ShellSet as a member, for easy retrieval of the shell indices.
 */
template <typename _Shell, size_t _N, typename _IntegralScalar>
class LibcintOneElectronIntegralEngine: public BaseOneElectronIntegralEngine<_Shell, _N, _IntegralScalar> {
public:
    using Shell = _Shell;                    // the type of shell the integral engine is able to handle
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;            // the number of components the operator has


private:
    Libcint1eFunction libcint_function;  // the libcint one-electron integral function

    // Data that has to be kept as a member (see the class note).
    libcint::RawContainer libcint_raw_container;  // the raw libcint data
    ShellSet<Shell> shell_set;                    // the corresponding shell set

    // Parameters to pass to the buffer.
    IntegralScalar scaling_factor = 1.0;  // a factor that is multiplied to all of the calculated integrals


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct an integral engine for the electronic dipole operator using libcint as a backend.
     * 
     *  @param op               the electronic electric dipole operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     */
    LibcintOneElectronIntegralEngine(const ElectronicDipoleOperator& op, const ShellSet<Shell>& shell_set) :
        libcint_function {LibcintInterfacer().oneElectronFunction(op)},
        libcint_raw_container {LibcintInterfacer().convert(shell_set)},
        shell_set {shell_set},
        scaling_factor {-1.0} {  // apply the minus sign which comes from the charge of the electrons -e

        LibcintInterfacer().setCommonOrigin(this->libcint_raw_container, op.reference());
    }

    /**
     *  Construct an integral engine for the kinetic operator using libcint as a backend.
     * 
     *  @param op               the kinetic operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     */
    LibcintOneElectronIntegralEngine(const KineticOperator& op, const ShellSet<Shell>& shell_set) :
        libcint_function {LibcintInterfacer().oneElectronFunction(op)},
        libcint_raw_container {LibcintInterfacer().convert(shell_set)},
        shell_set {shell_set} {}

    /**
     *  Construct an integral engine for the nuclear attraction operator using libcint as a backend.
     * 
     *  @param op               the nuclear attraction operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     */
    LibcintOneElectronIntegralEngine(const NuclearAttractionOperator& op, const ShellSet<Shell>& shell_set) :
        libcint_function {LibcintInterfacer().oneElectronFunction(op)},
        libcint_raw_container {LibcintInterfacer().convert(shell_set)},
        shell_set {shell_set} {}

    /**
     *  Construct an integral engine for the overlap operator using libcint as a backend.
     * 
     *  @param op               the overlap operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     */
    LibcintOneElectronIntegralEngine(const OverlapOperator& op, const ShellSet<Shell>& shell_set) :
        libcint_function {LibcintInterfacer().oneElectronFunction(op)},
        libcint_raw_container {LibcintInterfacer().convert(shell_set)},
        shell_set {shell_set} {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @param shell1           the first shell
     *  @param shell2           the second shell
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const GTOShell& shell1, const GTOShell& shell2) override {

        // Find to which indices in the RawContainer the given shells correspond
        int shell_indices[2];
        shell_indices[0] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell1));
        shell_indices[1] = static_cast<int>(findElementIndex(this->shell_set.asVector(), shell2));


        // Pre-allocate a raw buffer, because libcint functions expect a data pointer
        const size_t nbf1 = shell1.numberOfBasisFunctions();
        const size_t nbf2 = shell2.numberOfBasisFunctions();
        double libcint_buffer[2 * N * nbf1 * nbf2];


        // Let libcint compute the integrals and return the corresponding buffer
        const auto result = this->libcint_function(libcint_buffer, shell_indices, libcint_raw_container.atmData(), libcint_raw_container.numberOfAtoms(), libcint_raw_container.basData(), libcint_raw_container.numberOfBasisFunctions(), libcint_raw_container.envData());

        const std::vector<double> buffer_converted {libcint_buffer, libcint_buffer + N * nbf1 * nbf2};  // std::vector constructor from .begin() and .end()

        return std::make_shared<LibcintOneElectronIntegralBuffer<IntegralScalar, N>>(buffer_converted, nbf1, nbf2, result, this->scaling_factor);
    }
};


}  // namespace GQCP
