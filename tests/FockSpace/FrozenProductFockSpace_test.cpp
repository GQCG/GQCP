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
#define BOOST_TEST_MODULE "FrozenProductFockSpace"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "FockSpace/FrozenProductFockSpace.hpp"



BOOST_AUTO_TEST_CASE ( FrozenProductFockSpace_constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 1));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 2));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 3));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 4));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (10, 5, 5, 5));

    GQCP::ProductFockSpace product_fock_space (10, 5, 5);

    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 1));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 2));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 3));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 4));
    BOOST_CHECK_NO_THROW(GQCP::FrozenProductFockSpace (product_fock_space, 5));
}

BOOST_AUTO_TEST_CASE ( FrozenProductFockSpace_member_test ) {

    GQCP::FrozenProductFockSpace frozen_space (10, 5, 5, 2);

    const GQCP::FrozenFockSpace& alpha_member = frozen_space.get_frozen_fock_space_alpha();
    const GQCP::FrozenFockSpace& beta_member = frozen_space.get_frozen_fock_space_beta();

    BOOST_CHECK(alpha_member.get_N() == 5);
    BOOST_CHECK(beta_member.get_N() == 5);

    BOOST_CHECK(alpha_member.get_K() == 10);
    BOOST_CHECK(beta_member.get_K() == 10);
}
