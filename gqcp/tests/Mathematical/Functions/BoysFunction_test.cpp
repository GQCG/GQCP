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

#define BOOST_TEST_MODULE "BoysFunction"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/ChronusQ/engines.hpp"
#include "Mathematical/Functions/BoysFunction.hpp"


/**
 *  Check some values for the complex Boys function, with smaller values for z.
 * 
 *  The reference values are calculated through Wolfram Alpha.
 */
BOOST_AUTO_TEST_CASE(complex_boys_small) {

    const std::complex<double> z1 {1.0, 1.0};
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(0, z1) - std::complex<double>(0.69871273063423359570, -0.17867289388819850452)) < 1.0e-15);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(1, z1) - std::complex<double>(0.15770840051318358725, -0.09226490963072674028)) < 1.0e-15);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(2, z1) - std::complex<double>(0.07678105948851744967, -0.06039848610805146088)) < 1.0e-15);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(3, z1) - std::complex<double>(0.04817665805225730043, -0.04439293549582985343)) < 1.0e-15);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(4, z1) - std::complex<double>(0.03431995580042284671, -0.03491529220927123451)) < 1.0e-15);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(5, z1) - std::complex<double>(0.02635893440676594190, -0.02869781152193039798)) < 1.0e-15);
}


/**
 *  Check some values for the complex Boys function, with larger values for z.
 * 
 *  The reference values are calculated through Wolfram Alpha.
 */
BOOST_AUTO_TEST_CASE(complex_boys_large) {

    const std::complex<double> z2 {20.0, 12.0};
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(0, z2) - std::complex<double>(0.17684541398663813126, -0.04898334396254606203)) < 1.0e-12);
    std::cout << GQCP::BoysFunction()(0, z2) << std::endl;
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(1, z2) - std::complex<double>(0.00271057730158496915, -0.00285093000766361338)) < 1.0e-12);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(2, z2) - std::complex<double>(0.00005514807882485719, -0.00024690862551866566)) < 1.0e-12);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(3, z2) - std::complex<double>(-8.5475714213213795695e-06, -0.00002573506298602072)) < 1.0e-12);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(4, z2) - std::complex<double>(-3.0868137385792043864e-06, -2.65157542838644449488e-06)) < 1.0e-12);
    BOOST_CHECK(std::abs(GQCP::BoysFunction()(5, z2) - std::complex<double>(-7.7393812221822929702e-07, -1.32269247036353408729e-07)) < 1.0e-12);
}
