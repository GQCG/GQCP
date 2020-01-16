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
#define BOOST_TEST_MODULE "OrbitalRotationGenerators_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/OrbitalRotationGenerators.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    const GQCP::VectorX<double> kappa (4);  // not a triangular number
    BOOST_CHECK_THROW(GQCP::OrbitalRotationGenerators generators (kappa), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( asMatrix ) {

    GQCP::VectorX<double> kappa (3);
    kappa << 1, 2, 3;

    GQCP::SquareMatrix<double> kappa_matrix (3);
    kappa_matrix << 0, -1, -2,
                    1,  0, -3,
                    2,  3,  0;

    GQCP::OrbitalRotationGenerators generators (kappa);
    BOOST_CHECK(generators.asMatrix().isApprox(kappa_matrix));
}


BOOST_AUTO_TEST_CASE ( FromOccOcc ) {

    GQCP::VectorX<double> kappa (3);
    kappa << 1, 2, 3;
    GQCP::OrbitalRotationGenerators occ_occ_generators (kappa);  // 3 occupied spatial orbitals

    auto full_generators = GQCP::OrbitalRotationGenerators::FromOccOcc(occ_occ_generators, 4);  // 4 total spatial orbitals

    GQCP::SquareMatrix<double> full_kappa_matrix (4);
    full_kappa_matrix << 0, -1, -2, 0,
                         1,  0, -3, 0,
                         2,  3,  0, 0,
                         0,  0,  0, 0;
    BOOST_CHECK(full_generators.asMatrix().isApprox(full_kappa_matrix));
}
