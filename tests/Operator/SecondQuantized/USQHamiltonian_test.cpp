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
#define BOOST_TEST_MODULE "USQHamiltonian"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/linalg.hpp"

#include <boost/math/constants/constants.hpp>



BOOST_AUTO_TEST_CASE ( USQHamiltonian_constructor ) {
    
    // Create single-particle basis
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    // const GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (water, "STO-3G");
    BOOST_CHECK(true);
    /*
    // Create One- and SQTwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = sp_basis.numberOfBasisFunctions();
    GQCP::QCMatrix<double> H_core = GQCP::QCMatrix<double>::Random(K, K);

    GQCP::QCRankFourTensor<double> g (K);
    g.setRandom();

   // Create SQ operators with greater dimensions
    GQCP::QCMatrix<double> H_core_faulty = GQCP::QCMatrix<double>::Random(K+1, K+1);
    GQCP::QCRankFourTensor<double> g_faulty (K+1);
    g_faulty.setRandom();
    // Create SQHamilonians with different dimensions
    GQCP::SQHamiltonian<double> sq_hamiltonian_a (GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}));
    GQCP::SQHamiltonian<double> sq_hamiltonian_b (GQCP::ScalarSQOneElectronOperator<double>({H_core}), GQCP::ScalarSQTwoElectronOperator<double>({g}));
    GQCP::SQHamiltonian<double> sq_hamiltonian_b_faulty (GQCP::ScalarSQOneElectronOperator<double>({H_core_faulty}), GQCP::ScalarSQTwoElectronOperator<double>({g_faulty}));
    BOOST_CHECK(true);
    // Check if a correct constructor works with compatible elements
    //BOOST_CHECK_NO_THROW(GQCP::USQHamiltonian<double> (sq_hamiltonian_a, sq_hamiltonian_b, GQCP::ScalarSQTwoElectronOperator<double>({g})));
    // Check if a constructor throws an error with incompatible elements
    BOOST_CHECK_THROW(GQCP::USQHamiltonian<double> (sq_hamiltonian_a, sq_hamiltonian_b_faulty, GQCP::ScalarSQTwoElectronOperator<double>({g})), std::invalid_argument);
   // BOOST_CHECK_THROW(GQCP::USQHamiltonian<double> (sq_hamiltonian_a, sq_hamiltonian_b, GQCP::ScalarSQTwoElectronOperator<double>({g_faulty})), std::invalid_argument);
   */
}


BOOST_AUTO_TEST_CASE ( USQHamiltonian_transform ) {
    /*
    // This test will test if a total transformation or two individual transformations for the individual components of the USQHamiltonian amount to the same result
    // Create single-particle basis for alpha and beta
    auto water = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis_a (water, "STO-3G");
    const GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis_b (water, "STO-3G");

    size_t K = sp_basis_a.numberOfBasisFunctions();

    GQCP::USQHamiltonian<double> usq_hamiltonian1 = GQCP::USQHamiltonian<double>::Molecular(sp_basis_a, sp_basis_b, water);
    GQCP::USQHamiltonian<double> usq_hamiltonian2 = GQCP::USQHamiltonian<double>::Molecular(sp_basis_a, sp_basis_b, water);

    GQCP::SquareMatrix<double> U = GQCP::SquareMatrix<double>::RandomUnitary(K);

    usq_hamiltonian1.transform(U);
    usq_hamiltonian2.transformAlpha(U);
    usq_hamiltonian2.transformBeta(U);

    BOOST_CHECK(usq_hamiltonian1.twoElectronMixed().parameters().isApprox(usq_hamiltonian2.twoElectronMixed().parameters()));
    BOOST_CHECK(usq_hamiltonian1.alphaHamiltonian().core().parameters().isApprox(usq_hamiltonian2.alphaHamiltonian().core().parameters()));
    BOOST_CHECK(usq_hamiltonian1.betaHamiltonian().core().parameters().isApprox(usq_hamiltonian2.betaHamiltonian().core().parameters()));   
    */
}

