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


BOOST_AUTO_TEST_CASE(crawdad_sto3g_water) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt)
    double ref_energy_correction = -0.049149636120;


    // Create the molecular Hamiltonian in the RHF basis
    auto water = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {water, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, water);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(water.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.coefficientMatrix());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(sq_hamiltonian, rhf_parameters);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE(crawdad_sto3g_methane) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/ch4_sto3g/output.txt)
    double ref_energy_correction = -0.056046676165;

    // Create the molecular Hamiltonian in the RHF basis
    auto methane = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {methane, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, methane);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(methane.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    sq_hamiltonian.transform(rhf_parameters.coefficientMatrix());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(sq_hamiltonian, rhf_parameters);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}
