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
#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGPSESolver.hpp"
#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "Properties/expectation_values.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


/**
 *  Check if the analytical AP1roG energy is equal to the contraction of the 1- and 2-DM with the 1- and 2-electron integrals
 */
BOOST_AUTO_TEST_CASE ( energy_as_contraction ) {

    // Prepare the molecular Hamiltonian in the RHF basis
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const auto N_P = h2.numberOfElectrons()/2;
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, h2);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, sp_basis, h2);
    plain_scf_solver.solve();
    const auto rhf = plain_scf_solver.get_solution();
    GQCP::basisTransform(sp_basis, sq_hamiltonian, rhf.get_C());


    // Solve the AP1roG PSEs for the geminal coefficients
    GQCP::AP1roGPSEs pses (sq_hamiltonian, N_P);
    const GQCP::AP1roGPSESolver pse_solver (pses);
    const auto G = pse_solver.solve();  // initial guess is zero
    std::cout << "G: " << std::endl << G.asMatrix() << std::endl << std::endl;


    // Determine the Lagrangian multipliers in order to calculate the DMs
    GQCP::AP1roGLagrangianOptimizer lagrangian_optimizer (G, sq_hamiltonian);
    const auto multipliers = lagrangian_optimizer.solve();
    std::cout << "lambda: " << std::endl << multipliers.asMatrix() << std::endl << std::endl;


    // Calculate the DMs and check the trace with the one- and two-electron integrals
    const double electronic_energy = GQCP::calculateAP1roGEnergy(G, sq_hamiltonian);
    std::cout << "E: " << std::endl << electronic_energy << std::endl << std::endl;

    const auto D = GQCP::calculate1RDM(G, multipliers);
    const auto d = GQCP::calculate2RDM(G, multipliers);
    const double electronic_energy_by_contraction = GQCP::calculateExpectationValue(sq_hamiltonian, D, d);
    std::cout << "E_contr: " << std::endl << electronic_energy_by_contraction << std::endl << std::endl;

    BOOST_CHECK(std::abs(electronic_energy_by_contraction - electronic_energy) < 1.0e-09);
}
