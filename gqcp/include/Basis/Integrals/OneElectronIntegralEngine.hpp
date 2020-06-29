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
#include "Basis/Integrals/McMurchieDavidsonCoefficient.hpp"
#include "Basis/Integrals/OneElectronIntegralBuffer.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/CartesianDirection.hpp"

#include <boost/math/constants/constants.hpp>


namespace GQCP {


/**
 *  An integral engine that can calculate overlap integrals over shells.
 */
class OverlapIntegralEngine:
    public BaseOneElectronIntegralEngine<GTOShell, 1, double> {
public:
    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Calculate all the overlap integrals over the given shells.
     * 
     *  @param shell1           the first shell
     *  @param shell2           the second shell
     * 
     *  @note This method is not marked const to allow the Engine's internals to be changed
     * 
     *  @return a buffer containing the calculated integrals
     */
    std::shared_ptr<BaseOneElectronIntegralBuffer<IntegralScalar, N>> calculate(const Shell& shell1, const Shell& shell2) override {

        // Prepare some variables.
        const auto K = shell1.nucleus().position();
        const auto L = shell2.nucleus().position();

        const auto& gaussian_exponents1 = shell1.gaussianExponents();
        const auto& gaussian_exponents2 = shell2.gaussianExponents();

        const auto& contraction_coefficients1 = shell1.contractionCoefficients();
        const auto& contraction_coefficients2 = shell2.contractionCoefficients();


        // Loop over all basis functions that are contained in the shell. Since they are contracted GTOs, we will have to generate all the primitives.
        const auto all_cartesian_exponents1 = shell1.generateCartesianExponents();
        const auto all_cartesian_exponents2 = shell2.generateCartesianExponents();
        std::vector<double> integrals;  // a buffer that stores the calculated integrals
        for (const auto cartesian_exponents1 : all_cartesian_exponents1) {
            for (const auto cartesian_exponents2 : all_cartesian_exponents2) {

                // Calculate the contracted integral as a contraction over the contraction coefficients and the primitive integrals.
                double integral = 0.0;
                for (size_t c1 = 0; c1 < shell1.contractionSize(); c1++) {
                    const auto alpha = gaussian_exponents1[c1];
                    const auto d1 = contraction_coefficients1[c1];

                    for (size_t c2 = 0; c2 < shell2.contractionSize(); c2++) {
                        const auto beta = gaussian_exponents2[c2];
                        const auto d2 = contraction_coefficients2[c2];

                        // THIS IS THE INTEGRAL OVER TWO GAUSSIAN PRIMITIVES
                        // THIS BEHAVIOUR CAN BE OVERRIDDEN BY DIFFERENT INTEGRAL ENGINES
                        const auto p = alpha + beta;
                        double primitive_integral = std::pow(boost::math::constants::pi<double>() / p, 1.5);

                        for (const auto& direction : {GQCP::CartesianDirection::x, GQCP::CartesianDirection::y, GQCP::CartesianDirection::z}) {
                            const McMurchieDavidsonCoefficient E {K(direction), alpha, L(direction), beta};
                            const auto i = cartesian_exponents1.value(direction);
                            const auto j = cartesian_exponents2.value(direction);

                            primitive_integral *= E(i, j, 0);
                        }
                        // THIS IS THE INTEGRAL OVER TWO GAUSSIAN PRIMITIVES


                        integral += d1 * d2 * primitive_integral;
                    }
                }

                integrals.push_back(integral);
            }
        }

        std::array<std::vector<IntegralScalar>, N> buffer {integrals};
        return std::make_shared<OneElectronIntegralBuffer<double, 1>>(shell1.numberOfBasisFunctions(), shell2.numberOfBasisFunctions(), buffer);
    }
};


}  // namespace GQCP
