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
#define BOOST_TEST_MODULE "BasisFunction"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "Basis/BasisFunction.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero(3);
    GQCP::CartesianGTO gto1 (1.0, GQCP::CartesianExponents(1, 0, 0), center);
    GQCP::CartesianGTO gto2 (1.1, GQCP::CartesianExponents(1, 0, 0), center);
    GQCP::CartesianGTO gto3 (1.1, GQCP::CartesianExponents(0, 1, 0), center);

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 {{0.5, 1.0}, {gto1, gto2}};
    BOOST_CHECK_NO_THROW(GQCP::BasisFunction bf1 (lc1));

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 {{0.5, 1.0}, {gto1, gto3}};
    BOOST_CHECK_THROW(GQCP::BasisFunction bf2 (lc2), std::invalid_argument);
}



BOOST_AUTO_TEST_CASE ( operator_equals ) {

    GQCP::Vector<double, 3> center = GQCP::Vector<double, 3>::Zero(3);

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(2, 0, 0), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(2, 0, 0), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(1, 1, 0), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(1, 1, 0), center)}};


    GQCP::BasisFunction bf1 (lc1);
    GQCP::BasisFunction bf2 (lc2);

    BOOST_CHECK(bf1 == bf1);
    BOOST_CHECK(!(bf1 == bf2));
}
