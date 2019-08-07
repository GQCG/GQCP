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
#define BOOST_TEST_MODULE "OneRDM"


#include <boost/test/included/unit_test.hpp>

#include "RDM/OneRDM.hpp"


BOOST_AUTO_TEST_CASE ( OneRDM_constructor ) {

    // Check a correct constructor
    GQCP::MatrixX<double> matrix = GQCP::MatrixX<double>::Zero(4, 4);
    GQCP::OneRDM<double> D (matrix);


    // Check a faulty constructor
    GQCP::MatrixX<double> matrix2 = GQCP::MatrixX<double>::Zero(3, 4);
    BOOST_CHECK_THROW(GQCP::OneRDM<double> D2 (matrix2), std::invalid_argument);
}
