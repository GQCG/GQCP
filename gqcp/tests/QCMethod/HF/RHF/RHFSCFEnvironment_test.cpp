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

#define BOOST_TEST_MODULE "RHFSCFEnvironment"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"


/**
 *  Check if the RHFSCFEnvironment constructor throws when given an odd number of electrons.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    // Check a correct constructor with an even number of electrons
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {h2, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    const auto K = sq_hamiltonian.numberOfOrbitals();                                 // the number of spatial orbitals

    BOOST_CHECK_NO_THROW(const GQCP::RHFSCFEnvironment<double> rhf_environment(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters(),
                                                                               GQCP::RTransformation<double>::Random(K)));


    // Check if a faulty constructor with an odd number of electron throws
    auto h2_ion = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz", +1);
    BOOST_CHECK_THROW(const GQCP::RHFSCFEnvironment<double> rhf_environment(h2_ion.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters(),
                                                                            GQCP::RTransformation<double>::Random(K)),
                      std::invalid_argument);
}
