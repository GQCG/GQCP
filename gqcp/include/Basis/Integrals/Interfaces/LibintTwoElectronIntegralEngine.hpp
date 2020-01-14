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

#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralBuffer.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"


namespace GQCP {


/**
 *  A two-electron integral engine that uses libint as its backend
 * 
 *  @tparam _N                  the number of components the operator has
 */
template <size_t _N>
class LibintTwoElectronIntegralEngine : public BaseTwoElectronIntegralEngine<GTOShell, _N, double> {
public:
    using IntegralScalar = double;  // the scalar representation of an integral for libint is always a real number
    static constexpr auto N = _N;  // the number of components the operator has


private:
    libint2::Engine libint2_engine;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     */
    LibintTwoElectronIntegralEngine(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) :
        libint2_engine (LibintInterfacer::get().createEngine(op, max_nprim, max_l))
    {}



    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @param shell1           the first shell
     *  @param shell2           the second shell
     *  @param shell3           the third shell
     *  @param shell4           the fourth shell
     * 
     *  This method is not marked const to allow the Engine's internals to be changed
     */
    std::shared_ptr<BaseTwoElectronIntegralBuffer<IntegralScalar, N>> calculate(const GTOShell& shell1, const GTOShell& shell2, const GTOShell& shell3, const GTOShell& shell4) override {

        const auto libint_shell1 = LibintInterfacer::get().interface(shell1);
        const auto libint_shell2 = LibintInterfacer::get().interface(shell2);
        const auto libint_shell3 = LibintInterfacer::get().interface(shell3);
        const auto libint_shell4 = LibintInterfacer::get().interface(shell4);

        const auto& libint2_buffer = this->libint2_engine.results();
        this->libint2_engine.compute(libint_shell1, libint_shell2, libint_shell3, libint_shell4);
        return std::make_shared<LibintTwoElectronIntegralBuffer<N>>(libint2_buffer, shell1.numberOfBasisFunctions(), shell2.numberOfBasisFunctions(), shell3.numberOfBasisFunctions(), shell4.numberOfBasisFunctions());
    }
};


}  // namespace GQCP
