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
#define BOOST_TEST_MODULE "SQTwoElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"

#include "Utilities/linalg.hpp"
#include "Utilities/miscellaneous.hpp"


/**
 *  Check the interface for constructing SQTwoElectronOperators from Tensors
 */
BOOST_AUTO_TEST_CASE ( SQTwoElectronOperator_constructor ) {

    // Check a correct constructor
    const GQCP::QCRankFourTensor<double> tensor (3);
    GQCP::ScalarSQTwoElectronOperator<double> O {tensor};


    // Check a faulty constructor
    GQCP::Tensor<double, 4> tensor2 (3, 3, 3, 2);
    BOOST_CHECK_THROW(GQCP::ScalarSQTwoElectronOperator<double> O2 {tensor2}, std::invalid_argument);
}


/**
 *  Check if the zero constructor really sets is parameters to all zeros
 */
BOOST_AUTO_TEST_CASE ( SQTwoElectronOperator_zero_constructor ) {

    const size_t dim = 2;
    GQCP::ScalarSQOneElectronOperator<double> op {dim};

    BOOST_CHECK_EQUAL(op.dimension(), dim);
    BOOST_CHECK(op.parameters().isZero(1.0e-08));
}


/**
 *  Check if the formulas in effectiveOneElectronPartition are implemented correctly
 */
BOOST_AUTO_TEST_CASE ( SQTwoElectronOperator_effectiveOneElectronPartition ) {

    const size_t K = 4;
    auto K_ = static_cast<double>(K);

    // Set up toy 2-electron integrals
    GQCP::QCRankFourTensor<double> g_par (K);
    g_par.setZero();

    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < K; j++) {
            for (size_t k = 0; k < K; k++) {
                for (size_t l = 0; l < K; l++) {
                    g_par(i,j,k,l) = (i+1) + 2*(j+1) + 4*(k+1) + 8*(l+1);
                }
            }
        }
    }

    GQCP::ScalarSQTwoElectronOperator<double> g {g_par};


    // Set up the reference effective one-electron integrals by manual calculation
    GQCP::QCMatrix<double> k_par_ref = GQCP::QCMatrix<double>::Zero(K, K);  // reference parameters
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            auto p_ = static_cast<double>(p) + 1;
            auto q_ = static_cast<double>(q) + 1;

            k_par_ref(p,q) = -K_ / 2 * (p_ + 8*q_ + 3*K_ + 3);
        }
    }


    BOOST_CHECK(k_par_ref.isApprox(g.effectiveOneElectronPartition().parameters(), 1.0e-08));
}
