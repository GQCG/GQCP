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
#define BOOST_TEST_MODULE "WaveFunction"

#include <boost/test/unit_test.hpp>

#include "FockSpace/WaveFunction/SelectedWaveFunction.hpp"

#include "Basis/transform.hpp"
#include "FockSpace/FockSpace.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


/*
 *  Tests the non-trivial constructor of the SelectedWaveFunction by selecting the 2 largest coefficients
 */ 
BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCP::ONV onv1 = GQCP::ONV(3, 1, 1);  // 001
    GQCP::ONV onv2 = GQCP::ONV(3, 1, 2);  // 010

    Eigen::Vector4d test_coefficients;
    test_coefficients << 2.0/7, -5.0/7, 4.0/7, 2.0/7;

    GQCP::SelectedFockSpace selected_fock_space (3, 1, 1);

    selected_fock_space.addConfiguration(onv1.asString(), onv1.asString());  // 001 | 001
    selected_fock_space.addConfiguration(onv1.asString(), onv2.asString());  // 001 | 010
    selected_fock_space.addConfiguration(onv2.asString(), onv1.asString());  // 010 | 001
    selected_fock_space.addConfiguration(onv2.asString(), onv2.asString());  // 010 | 010

    // Create wavefunction
    GQCP::WaveFunction wf (selected_fock_space, test_coefficients);

    // Create selected wavefunction with 2 of the largest coefficients (configuration 2 and 3)
    GQCP::SelectedWaveFunction wf2 (wf, 2);
    const auto& configurations = wf2.fockSpace().get_configurations();
    const auto& coefficients = wf2.get_coefficients();
    
    // Given the extraction method the smallest coefficient will come first this is configuration 3: 010 | 001
    BOOST_CHECK(configurations[0].onv_alpha.asString() == "010");
    BOOST_CHECK(configurations[0].onv_beta.asString() == "001");
    BOOST_CHECK(std::abs(coefficients(0) - (4.0/7)) < 1.0e-9);

    // Followed by configuration 2: 001 | 010
    BOOST_CHECK(configurations[1].onv_alpha.asString() == "001");
    BOOST_CHECK(configurations[1].onv_beta.asString() == "010");
    BOOST_CHECK(std::abs(coefficients(1) - (-5.0/7)) < 1.0e-9);
}
