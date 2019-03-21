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

    std::vector<double> exp {1.0, 1.1};
    std::vector<double> coeff {0.5, 1.0};

    GQCP::Shell s_shell (0, GQCP::Atom(), exp, coeff);
    GQCP::Shell p_shell (1, GQCP::Atom(), exp, coeff);
    GQCP::Shell d_shell (2, GQCP::Atom(), exp, coeff);


//    for (const auto& bf : s_shell.basisFunctions()) {
//        std::cout << "Basis function: " << std::endl;
//        for (const auto& func : bf.get_functions()) {
//            std::cout << "\tGaussian exponent: " << func.get_gaussian_exponent() << std::endl;
//        }
//    }

    for (const auto& bf : p_shell.basisFunctions()) {
        std::cout << "Basis function: " << std::endl;
        for (const auto& func : bf.get_functions()) {
            std::cout << "\tGaussian exponent: " << func.get_gaussian_exponent() << std::endl;
        }
    }

}
