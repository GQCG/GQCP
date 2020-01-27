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
#define BOOST_TEST_MODULE "ONV"

#include <boost/test/unit_test.hpp>

#include "FockSpace/Configuration.hpp"


BOOST_AUTO_TEST_CASE ( CompactConfigurationString ) {

    // Create ONV : 10 considered bits and 5 set bits with distributed as "0000011111" = 31
    GQCP::ONV onv1 (10, 5, 31);  // "0000011111" = 31
    GQCP::ONV onv2 (10, 5, 62);  // "0000111110" = 62

    std::string compact_reference = "0000122221";
    std::string compact_reference_delimiter = "0 0 0 0 1 2 2 2 2 1";

    GQCP::Configuration configuration {onv1, onv2};

    std::string compact_str = configuration.spinSummedRepresentation();
    std::string compact_str_delimiter = configuration.spinSummedRepresentation(" ");

    BOOST_CHECK(compact_str.compare(compact_reference) == 0);
    BOOST_CHECK(compact_str_delimiter.compare(compact_reference_delimiter) == 0);
}

