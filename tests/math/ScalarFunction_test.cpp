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
#define BOOST_TEST_MODULE "ScalarFunction"

#include <boost/test/unit_test.hpp>
//#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "math/ScalarFunction.hpp"

#include "Basis/CartesianGTO.hpp"



BOOST_AUTO_TEST_CASE ( ScalarFunctionProduct ) {

    // LinearCombination inherits from ScalarFunction, so we can use the class to test the behavior of ScalarFunctionProduct
    // Therefore, this test also makes sure that ScalarFunctionProduct correctly compiles with T1 and T2 having the same ::Cols and ::Scalar

    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero();

    GQCP::CartesianGTO gto1 (1.0, {1, 0, 0}, center);
    double coefficient1 = 1.0;
    double N1 = gto1.calculateNormalizationFactor();

    GQCP::CartesianGTO gto2 (2.00, {2, 0, 0}, center);
    double coefficient2 = -2.0;
    double N2 = gto2.calculateNormalizationFactor();

    GQCP::CartesianGTO gto3 (0.50, {0, 0, 1}, center);
    double coefficient3 = 0.75;
    double N3 = gto3.calculateNormalizationFactor();


    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 (coefficient1, gto1);
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 ({coefficient2, coefficient3}, {gto2, gto3});


    // Evaluate the product of two LinearCombinations at a point
    GQCP::Vector<double, 3> r = GQCP::Vector<double, 3>::Zero();
    r << 1.0, 1.0, 1.0;

    double ref_evaluation = 1*N1*std::exp(-3.0) * (-2*N2*std::exp(-6.0) + 0.75*N3*std::exp(-1.5));  // manual calculation
    double evaluation = (lc1 * lc2)(r);


    BOOST_CHECK(std::abs(evaluation - ref_evaluation) < 1.0e-12);
}
