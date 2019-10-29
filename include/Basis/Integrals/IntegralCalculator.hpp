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


#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"


#include <array>
#include <memory>


namespace GQCP {


/**
 *  A class that calculates integrals over ShellSets: it loops over all shells in the given ShellSets. It just has two static member functions, for cleaner interface calls.
 */
class IntegralCalculator {
public:

    /*
     *  PUBLIC METHODS
     */

    /**
     *  Calculate all one-electron integrals over the basis functions inside the given ShellSets
     * 
     *  @tparam Shell                   the type of shell the integral engine is able to handle
     *  @tparam N                       the number of components the operator has
     *  @tparam IntegralScalar          the scalar representation of an integral
     *  @param engine                   the engine that can calculate one-electron integrals over shells (not const because we allow for non-const Engine::calculate() calls)
     *  @param shell_set                the set of shells over which the integrals should be calculated
     */
    template <typename Shell, size_t N, typename IntegralScalar>
    static auto calculate(BaseOneElectronIntegralEngine<Shell, N, IntegralScalar>& engine, const ShellSet<Shell>& shell_set) -> std::array<QCMatrix<IntegralScalar>, N> {

        // Initialize the N components of the matrix representations of the operator
        const auto nbf = shell_set.numberOfBasisFunctions();
        std::array<QCMatrix<IntegralScalar>, N> components;
        for (auto& component : components) {
            component = QCMatrix<IntegralScalar>::Zero(nbf, nbf);
        }


        // Loop over all shells (twice) inside the Shell set and let the engine calculate the integrals of two shells
        const auto nsh = shell_set.numberOfShells();
        const auto shells = shell_set.asVector();
        for (size_t sh1_index = 0; sh1_index < nsh; sh1_index++) {  // shell 1
            const auto bf1 = shell_set.basisFunctionIndex(sh1_index);
            const auto shell1 = shells[sh1_index];

            for (size_t sh2_index = 0; sh2_index < nsh; sh2_index++) {  // shell 2
                const auto bf2 = shell_set.basisFunctionIndex(sh2_index);
                const auto shell2 = shells[sh2_index];

                const auto buffer = engine.calculate(shell1, shell2);

                // Only if the integrals are not all zero, place them inside the full matrices
                if (buffer->areIntegralsAllZero()) {
                    continue;
                }
                buffer->emplace(components, bf1, bf2);
            }  // sh2_index
        }  // sh1_index

        return components;
    }


    /**
     *  Calculate all two-electron integrals over the basis functions inside the given ShellSets
     * 
     *  @tparam Shell                   the type of shell the integral engine is able to handle
     *  @tparam N                           the number of components the operator has
     *  @tparam IntegralScalar              the scalar representation of an integral
     * 
     *  @param engine                       the engine that can calculate two-electron integrals over shells
     *  @param shell_set                    the set of shells over which the integrals should be calculated
     */
    template <typename Shell, size_t N, typename IntegralScalar>
    static auto calculate(BaseTwoElectronIntegralEngine<Shell, N, IntegralScalar>& engine, const ShellSet<Shell>& shell_set) -> std::array<QCRankFourTensor<IntegralScalar>, N> {

        // Initialize the N components of the matrix representations of the operator
        const auto nbf = shell_set.numberOfBasisFunctions();
        std::array<QCRankFourTensor<IntegralScalar>, N> components;
        for (auto& component : components) {
            component = QCRankFourTensor<IntegralScalar>(nbf);
            component.setZero();
        }


        // Loop over all shells (four times) inside the Shell set and let the engine calculate the integrals over four shells
        const auto nsh = shell_set.numberOfShells();
        const auto shells = shell_set.asVector();
        for (size_t sh1_index = 0; sh1_index < nsh; sh1_index++) {  // shell 1
            const auto bf1 = shell_set.basisFunctionIndex(sh1_index);
            const auto shell1 = shells[sh1_index];

            for (size_t sh2_index = 0; sh2_index < nsh; sh2_index++) {  // shell 2
                const auto bf2 = shell_set.basisFunctionIndex(sh2_index);
                const auto shell2 = shells[sh2_index];

                for (size_t sh3_index = 0; sh3_index < nsh; sh3_index++) {  // shell 3
                    const auto bf3 = shell_set.basisFunctionIndex(sh3_index);
                    const auto shell3 = shells[sh3_index];

                    for (size_t sh4_index = 0; sh4_index < nsh; sh4_index++) {  // shell 4
                        const auto bf4 = shell_set.basisFunctionIndex(sh4_index);
                        const auto shell4 = shells[sh4_index];

                        const auto buffer = engine.calculate(shell1, shell2, shell3, shell4);

                        // Only if the integrals are not all zero, place them inside the full matrices
                        if (buffer->areIntegralsAllZero()) {
                            continue;
                        }
                        buffer->emplace(components, bf1, bf2, bf3, bf4);  // place the calculated integrals inside the full tensors
                    }  // sh4_index
                }  // sh3_index
            }  // sh2_index
        }  // sh1_index

        return components;
    }
};


}  // namespace GQCP
