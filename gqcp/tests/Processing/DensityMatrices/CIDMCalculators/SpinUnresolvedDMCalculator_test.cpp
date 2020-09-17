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

#define BOOST_TEST_MODULE "UnresolvedCIDMCalculator_test"

#include <boost/test/unit_test.hpp>

#include "Processing/DensityMatrices/CIDMCalculators/SpinUnresolvedDMCalculator.hpp"


BOOST_AUTO_TEST_CASE(throw_1and2_DMs) {

    // Create a test wave function
    size_t M = 5;
    size_t N = 4;
    GQCP::SpinUnresolvedONVBasis fock_space {M, N};

    GQCP::VectorX<double> coeff {fock_space.dimension()};
    coeff << 1, 1, -2, 4, -5;


    // not implemented yet and should throw
    GQCP::SpinUnresolvedDMCalculator d(fock_space);
    BOOST_CHECK_THROW(d.calculate1DM(coeff), std::runtime_error);
    BOOST_CHECK_THROW(d.calculate2DM(coeff), std::runtime_error);
}
