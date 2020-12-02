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

#define BOOST_TEST_MODULE "Simple1DM"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Test if the expectation value of a one-electron operator in different orbital bases is the same.
 */
BOOST_AUTO_TEST_CASE(one_electron_operator_expectation_value_different_orbital_bases) {

    // Prepare the molecular Hamiltonian in the AO basis, in order to proceed with an RHF SCF calculation.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);  // In the AO basis.
    const auto K = hamiltonian.numberOfOrbitals();

    // Do the RHF SCF calculation to retrieve the RHF MOs.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), hamiltonian, S.parameters());
    auto diis_rhf_scf_solver = GQCP::RHFSCFSolver<double>::DIIS();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {hamiltonian};

    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, diis_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Prepare three one-electron operators in different orbital basis.
    const auto& h_core_AO = hamiltonian.core();
    const auto h_core_MO = h_core_AO.transformed(rhf_parameters.expansion());

    const auto T_random = GQCP::RTransformation<double>::Random(K);
    const auto h_core_random = h_core_AO.transformed(T_random);

    // Prepare three density matrices in the corresponding orbital bases.
    const auto D_AO = rhf_parameters.calculateOrthonormalBasis1DM();
    const auto D_MO = D_AO.transformed(rhf_parameters.expansion());
    const auto D_random = D_AO.transformed(T_random);


    // Check if the expectation values match.
    const double exp_val_AO = h_core_AO.calculateExpectationValue(D_AO);
    const double exp_val_MO = h_core_MO.calculateExpectationValue(D_MO);
    const double exp_val_random = h_core_random.calculateExpectationValue(D_random);

    BOOST_CHECK(std::abs(exp_val_AO - exp_val_MO) < 1.0e-12);
    BOOST_CHECK(std::abs(exp_val_AO - exp_val_random) < 1.0e-12);
}
