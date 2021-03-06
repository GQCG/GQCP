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

#define BOOST_TEST_MODULE "IterativeIdentitiesHessianModifier_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"

#include <Eigen/Dense>


/**
 *  Check if an iterative identities Hessian modifier produces a positive definite Hessian.
 */
BOOST_AUTO_TEST_CASE(becomes_positive_definite) {

    // Provide an indefinite matrix. Its eigenvalues are approximately. -1.7 and 4.7
    GQCP::SquareMatrix<double> A {2};
    // clang-format off
    A << -1.0, -2.0,
         -2.0,  4.0;
    // clang-format on


    // Create a Hessian modifier and supply the indefinite matrix.
    GQCP::IterativeIdentitiesHessianModifier hessian_modifier {};  // Use default values for the parameters.
    const auto modified_hessian = hessian_modifier(A);

    // Check if the modified Hessian is positive definite.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver {modified_hessian};
    BOOST_CHECK(solver.eigenvalues().minCoeff() > 0);  // A positive definite Hessian has all positive eigenvalues!
}
