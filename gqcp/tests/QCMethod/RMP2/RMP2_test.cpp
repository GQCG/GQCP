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

#define BOOST_TEST_MODULE "RMP2"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCMethod/RMP2/RMP2.hpp"


/**
 *  Check the RMP2 energy correction with results from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt).
 *  The test system is H2O in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(crawdad_sto3g_H2O) {

    const double ref_energy_correction = -0.049149636120;


    // Create the molecular Hamiltonian in the RHF basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    hamiltonian.transform(rhf_parameters.expansion());  // Now in the RHF orbital basis.


    // Check if the RMP2 energy correction is correct.
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(hamiltonian, rhf_parameters);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}


/**
 *  Check the RMP2 energy correction with results from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt).
 *  The test system is CH4 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(crawdad_sto3g_CH4) {

    const double ref_energy_correction = -0.056046676165;

    // Create the molecular Hamiltonian in the RHF basis.
    auto molecule = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, spin_orbital_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    hamiltonian.transform(rhf_parameters.expansion());


    // Check if the RMP2 correction is correct.
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(hamiltonian, rhf_parameters);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}
