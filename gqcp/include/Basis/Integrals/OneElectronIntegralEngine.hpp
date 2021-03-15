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
#include "Basis/Integrals/OneElectronIntegralBuffer.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"


namespace GQCP {


/**
 *  An integral engine that can calculate one-electron integrals over shells.
 * 
 *  @tparam _PrimitiveIntegralEngine            The type of integral engine that is used for calculating integrals over primitives.
 */
template <typename _PrimitiveIntegralEngine>
class OneElectronIntegralEngine:
    public BaseOneElectronIntegralEngine<typename _PrimitiveIntegralEngine::Shell, _PrimitiveIntegralEngine::Components, typename _PrimitiveIntegralEngine::IntegralScalar> {
public:
    // The type of integral engine that is used for calculating integrals over primitives.
    using PrimitiveIntegralEngine = _PrimitiveIntegralEngine;

    // The type of shell that this engine can calculate integrals over.
    using Shell = typename _PrimitiveIntegralEngine::Shell;

    // The scalar representation of one of the integrals.
    using IntegralScalar = typename _PrimitiveIntegralEngine::IntegralScalar;

    // The number of components the operator has.
    static constexpr auto N = _PrimitiveIntegralEngine::Components;


private:
    // The integral engine that is used for calculating integrals over primitives.
    PrimitiveIntegralEngine primitive_engine;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param primitive_engine             The integral engine that is used for calculating integrals over primitives.
     */
    OneElectronIntegralEngine(const PrimitiveIntegralEngine& primitive_engine) :
        primitive_engine {primitive_engine} {}


    /*
     *  MARK: Integral calculations.
     */

    /**
     *  Calculate all the integrals over the given shells.
     * 
     *  @param shell1           The first shell.
     *  @param shell2           The second shell.
     * 
     *  @note This method is not marked const to allow the Engine's internals to be changed.
     * 
     *  @return A buffer containing the calculated integrals.
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const Shell& shell1, const Shell& shell2) override {

        // In this function, we loop over all basis functions that the shells contain.
        const auto basis_functions1 = shell1.basisFunctions();
        const auto basis_functions2 = shell2.basisFunctions();

        std::array<std::vector<IntegralScalar>, N> integrals;  // A "buffer" that stores the calculated integrals.

        for (size_t i = 0; i < N; i++) {  // Loop over all components of the operator.
            this->primitive_engine.prepareStateForComponent(i);

            for (const auto& bf1 : basis_functions1) {
                const auto& coefficients1 = bf1.coefficients();
                const auto& primitives1 = bf1.functions();

                for (const auto& bf2 : basis_functions2) {
                    const auto& coefficients2 = bf2.coefficients();
                    const auto& primitives2 = bf2.functions();

                    IntegralScalar integral {};
                    for (size_t c1 = 0; c1 < bf1.length(); c1++) {
                        const auto& d1 = coefficients1[c1];
                        const auto& primitive1 = primitives1[c1];

                        for (size_t c2 = 0; c2 < bf2.length(); c2++) {
                            const auto& d2 = coefficients2[c2];
                            const auto& primitive2 = primitives2[c2];

                            const auto primitive_integral = this->primitive_engine.calculate(primitive1, primitive2);
                            integral += d1 * d2 * primitive_integral;
                        }
                    }
                    integrals[i].push_back(integral);
                }
            }
        }

        return std::make_shared<OneElectronIntegralBuffer<IntegralScalar, N>>(shell1.numberOfBasisFunctions(), shell2.numberOfBasisFunctions(), integrals);
    }
};


}  // namespace GQCP
