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
#define BOOST_TEST_MODULE "AP1roG"

#include <boost/test/unit_test.hpp>

#include "Basis/transform.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"
#include "QCMethod/Geminals/PSEnvironment.hpp"
#include "QCMethod/Geminals/AP1roGLagrangianOptimizer.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "QCModel/Geminals/AP1roG.hpp"


/**
 *  Check if the analytical AP1roG energy is equal to the contraction of the 1- and 2-DM with the 1- and 2-electron integrals.
 */
BOOST_AUTO_TEST_CASE ( energy_as_contraction ) {

    // Prepare the molecular Hamiltonian in the AO basis.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const auto N_P = h2.numberOfElectrons()/2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis

    // Transform the Hamiltonian to the canonical RHF basis.
    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();
    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Optimize the AP1roG wave function model for the geminal coefficients.
    auto solver = GQCP::NonLinearEquationSolver<double>::Newton();
    auto environment = GQCP::PSEnvironment::AP1roG(sq_hamiltonian, N_P);  // zero initial guess
    const auto qc_structure = GQCP::QCMethod::AP1roG(sq_hamiltonian, N_P).optimize(solver, environment);

    const auto G = qc_structure.groundStateParameters().geminalCoefficients();
    const auto electronic_energy = qc_structure.groundStateEnergy();


    // Determine the Lagrangian multipliers in order to calculate the DMs.
    GQCP::AP1roGLagrangianOptimizer lagrangian_optimizer (G, sq_hamiltonian);
    const auto multipliers = lagrangian_optimizer.solve();


    // Calculate the DMs and check the trace with the one- and two-electron integrals.
    const auto D = GQCP::QCModel::AP1roG::calculate1RDM(G, multipliers);
    const auto d = GQCP::QCModel::AP1roG::calculate2RDM(G, multipliers);
    const double electronic_energy_by_contraction = sq_hamiltonian.calculateExpectationValue(D, d);

    BOOST_CHECK(std::abs(electronic_energy_by_contraction - electronic_energy) < 1.0e-09);
}
