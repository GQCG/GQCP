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

#define BOOST_TEST_MODULE "Function"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/Function.hpp"


/**
 *  `EvaluableLinearCombination` inherits from `Function`, so we can use the class to test the behavior of FunctionProduct.
 *  Therefore, this test also makes sure that `FunctionProduct` correctly compiles with T1 and T2 having the same InputType.
 */
BOOST_AUTO_TEST_CASE(FunctionProduct) {

    // Set up a linear combination of 3 CartesianGTOs
    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    const GQCP::CartesianGTO gto1 {1.0, {1, 0, 0}, center};
    const double coefficient1 = 1.0;

    const GQCP::CartesianGTO gto2 {2.00, {2, 0, 0}, center};
    const double coefficient2 = -2.0;

    const GQCP::CartesianGTO gto3 {0.50, {0, 0, 1}, center};
    const double coefficient3 = 0.75;


    const GQCP::EvaluableLinearCombination<double, GQCP::CartesianGTO> lc1 {coefficient1, gto1};
    const GQCP::EvaluableLinearCombination<double, GQCP::CartesianGTO> lc2 {{coefficient2, coefficient3}, {gto2, gto3}};


    // Evaluate the product of two EvaluableLinearCombinations at a position in space.
    GQCP::Vector<double, 3> r {1.0, 1.0, 1.0};
    BOOST_CHECK(std::abs((lc1 * lc2)(r) -0.0080849277955083707311) < 1.0e-12);  // reference value by WolframAlpha
}
