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
#ifndef GQCP_INTEGRALCALCULATOR_HPP
#define GQCP_INTEGRALCALCULATOR_HPP


#include "Basis/ShellSet.hpp"
#include "Basis/BaseOneElectronIntegralEngine.hpp"
#include "Mathematical/SquareMatrix.hpp"

#include <array>
#include <memory>


namespace GQCP {


/**
 *  A class that calculates integrals over ShellSets: it loops over all shells in the given ShellSets
 */
class IntegralCalculator {
public:

    /**
     *  Calculate all integrals over the basis functions inside the given ShellSets
     * 
     *  @param engine               the engine that can calculate integrals over shells
     *  @param shell_set            the set of shells over which the integrals should be calculated
     * 
     *  @tparam ShellType           the type of shell the integral engine is able to handle
     *  @tparam N                   the number of components the operator has
     *  @tparam Scalar              the scalar representation of an integral
     */
    template <typename ShellType, size_t N, typename Scalar>
    static auto calculate(const BaseOneElectronIntegralEngine<ShellType, N, Scalar>& engine, const ShellSet<ShellType>& shell_set) -> std::array<SquareMatrix<Scalar>, N> {

        // Initialize the N components of the matrix representations of the operator
        const auto nbf = shell_set.numberOfBasisFunctions();
        std::array<SquareMatrix<Scalar>, N> components;
        for (auto& component : components) {
            component = SquareMatrix<Scalar>::Zero(nbf, nbf);
        }


        // Loop over all shells inside the Shell set and let the engine calculate the integrals
        const auto nsh = shell_set.numberOfShells();
        for (size_t sh1 = 0; sh1 < nsh; sh2++) {  // shell 1
            const auto bf1 = shell_set.basisFunctionIndex(sh1);

            for (size_t sh2 = 0; sh2 < nsh; sh2++) {  // shell 2
                const auto bf2 = shell_set.basisFunctionIndex(sh2);

                const auto buffer = engine.calculate(sh1, sh2);  // calculate the integrals over the two shells
                buffer.emplace(components, bf1, bf2);  // place the calculated integrals inside the full matrices
            }  // sh1
        }  // sh2

        return components;
    }

    
};


}  // namespace GQCP



#endif  // GQCP_INTEGRALCALCULATOR_HPP
