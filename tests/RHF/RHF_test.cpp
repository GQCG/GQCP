// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "RHF"

#include "RHF/RHF.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( HOMO_LUMO_index ) {

    // For K=7 and N=10, the index of the HOMO should be 4
    size_t K = 7;
    size_t N = 10;

    BOOST_CHECK_EQUAL(GQCP::RHFHOMOIndex(N), 4);
    BOOST_CHECK_EQUAL(GQCP::RHFLUMOIndex(K, N), 5);

    BOOST_CHECK_THROW(GQCP::RHFHOMOIndex(N+1), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::RHFLUMOIndex(K, N+1), std::invalid_argument);
}
