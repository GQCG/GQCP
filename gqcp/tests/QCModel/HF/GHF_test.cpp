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

#define BOOST_TEST_MODULE "QCModel::GHF"

#include <boost/test/unit_test.hpp>

#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCModel/HF/GHF.hpp"


/**
 *  Check if the GHF energy is equal to the expectation value of the Hamiltonian through its density matrices.
 */
BOOST_AUTO_TEST_CASE(GHF_DMs) {

    // Perform a GHF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    auto hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(spinor_basis, molecule);  // In an AO basis.

    auto ghf_environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spinor_basis.overlap().parameters());
    auto plain_ghf_scf_solver = GQCP::GHFSCFSolver<double>::Plain();

    const auto ghf_qc_structure = GQCP::QCMethod::GHF<double>().optimize(plain_ghf_scf_solver, ghf_environment);
    const auto ghf_parameters = ghf_qc_structure.groundStateParameters();
    const auto ghf_energy = ghf_qc_structure.groundStateEnergy();

    // Determine the ghf energy through the expectation value of the Hamiltonian, and check the result.
    // Do the calculations in the GHF MO basis, in order to check the implementation of the GHF density matrices in MO basis.
    hamiltonian.transform(ghf_parameters.expansion());
    const auto D_MO = ghf_parameters.calculateOrthonormalBasis1DM();
    const auto d_MO = ghf_parameters.calculateOrthonormalBasis2DM();
    const double expectation_value = hamiltonian.calculateExpectationValue(D_MO, d_MO);

    BOOST_CHECK(std::abs(ghf_energy - expectation_value) < 1.0e-12);
}
