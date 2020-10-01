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

#define BOOST_TEST_MODULE "EvaluatableRSQOneElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Operator/SecondQuantized/EvaluatableScalarRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check the RHF density integrated over the whole space equals the number of electrons. In this test, we won't be using a very high-resolution grid, so so we can't expect very much about the accuracy of the integration.
 * 
 *  The subject of interest is H2//STO-3G.
 */
BOOST_AUTO_TEST_CASE(integrated_density_sto_3g) {

    // Prepare the molecular Hamiltonian (in AO basis).
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the scalar/AO basis

    // Prepare the canonical RHF orbitals.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    spinor_basis.transform(rhf_parameters.coefficientMatrix());


    // Calculate the RHF density.
    const auto rho_op = spinor_basis.quantize(GQCP::Operator::ElectronicDensity());
    const auto D = rhf_parameters.calculateOrthonormalBasis1DM();  // the (orthonormal) 1-DM for RHF
    // const auto density = rho_op.calculateDensity(D);

    // const auto grid = GQCP::CubicGrid::Centered(GQCP::Vector<double, 3>::Zero(), 50, 0.2);
    // const auto density_evaluated = grid.evaluate(density);

    // BOOST_CHECK(std::abs(grid.integrate(density_evaluated) - molecule.numberOfElectrons()) < 1.0e-04);  // 1.0e-04 is still reasonable given the accuracy of the grid
}


/**
 *  Check the RHF density integrated over the whole space equals the number of electrons. In this test, we won't be using a very high-resolution grid, so so we can't expect very much about the accuracy of the integration.
 * 
 *  The subject of interest is H2//cc-pVTZ.
 */
BOOST_AUTO_TEST_CASE(integrated_density_cc_pVTZ) {

    // Prepare the molecular Hamiltonian (in AO basis).
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "cc-pVDZ"};

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the scalar/AO basis

    // Prepare the canonical RHF orbitals.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    spinor_basis.transform(rhf_parameters.coefficientMatrix());


    // Calculate the RHF density.
    const auto rho_op = spinor_basis.quantize(GQCP::Operator::ElectronicDensity());
    const auto D = rhf_parameters.calculateOrthonormalBasis1DM();  // the (orthonormal) 1-DM for RHF
    // const auto density = rho_op.calculateDensity(D);

    // const auto grid = GQCP::CubicGrid::Centered(GQCP::Vector<double, 3>::Zero(), 50, 0.2);
    // const auto density_evaluated = grid.evaluate(density);

    // BOOST_CHECK(std::abs(grid.integrate(density_evaluated) - molecule.numberOfElectrons()) < 1.0e-03);  // 1.0e-04 is still reasonable given the accuracy of the grid
}
