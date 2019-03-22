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
#define BOOST_TEST_MODULE "Shell"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "Basis/Shell.hpp"


BOOST_AUTO_TEST_CASE ( Shell_constructor_throws ) {

    std::vector<double> exp1 {1.0, 1.1};
    std::vector<double> coeff1 {0.5, 1.0};
    std::vector<double> coeff2 {0.5, 1.0, 1.5};

    GQCP::Shell shell1 (0, GQCP::Atom(), exp1, coeff1);
    BOOST_CHECK_THROW(GQCP::Shell shell2 (0, GQCP::Atom(), exp1, coeff2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( basisFunctions ) {

    // Create a d-shell consisting of 2 contractions and test its basis functions
    std::vector<double> exp {1.0, 1.1};
    std::vector<double> coeff {0.5, 1.0};
    GQCP::Atom atom {};
    GQCP::Shell d_shell (2, atom, exp, coeff);

    auto center = atom.position;
    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc1 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(2, 0, 0), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(2, 0, 0), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc2 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(1, 1, 0), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(1, 1, 0), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc3 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(1, 0, 1), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(1, 0, 1), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc4 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(0, 2, 0), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(0, 2, 0), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc5 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(0, 1, 1), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(0, 1, 1), center)}};

    GQCP::LinearCombination<double, GQCP::CartesianGTO> lc6 {{0.5, 1.0},
        {GQCP::CartesianGTO(1.0, GQCP::CartesianExponents(0, 0, 2), center),
         GQCP::CartesianGTO(1.1, GQCP::CartesianExponents(0, 0, 2), center)}};


    GQCP::BasisFunction bf1 (lc1);
    GQCP::BasisFunction bf2 (lc2);
    GQCP::BasisFunction bf3 (lc3);
    GQCP::BasisFunction bf4 (lc4);
    GQCP::BasisFunction bf5 (lc5);
    GQCP::BasisFunction bf6 (lc6);

    std::vector<GQCP::BasisFunction> ref_bfs {bf1, bf2, bf3, bf4, bf5, bf6};

    BOOST_CHECK(ref_bfs == d_shell.basisFunctions());
}
