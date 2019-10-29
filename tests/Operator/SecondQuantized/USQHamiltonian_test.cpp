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




BOOST_AUTO_TEST_CASE ( USQHamiltonian_transform ) {

    // This test will test if a total transform or two individual transformations for the individual components of the unrestricted SQHamiltonian amount to the same result
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

    BOOST_CHECK(GQCP::SQHamiltonian<double> (usq_hamiltonian1.twoElectronMixed().parameters().isApprox(usq_hamiltonian2.twoElectronMixed().parameters()));
    BOOST_CHECK(GQCP::SQHamiltonian<double> (usq_hamiltonian1.alphaHamiltonian().core().parameters().isApprox(usq_hamiltonian2.alphaHamiltonian().core().parameters()));
    BOOST_CHECK(GQCP::SQHamiltonian<double> (usq_hamiltonian1.betaHamiltonian().core().parameters().isApprox(usq_hamiltonian2.betaHamiltonian().core().parameters()));

}

