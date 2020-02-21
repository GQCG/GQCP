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
#define BOOST_TEST_MODULE "AP1roGPSESolver"

#include <boost/test/unit_test.hpp>

#include "QCMethod/Geminals/AP1roGPSESolver.hpp"

#include "Basis/transform.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/Geminals/AP1roGPSEs.hpp"
#include "QCMethod/HF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF.hpp"
#include "QCMethod/HF/RHFSCFSolver.hpp"
#include "QCModel/Geminals/AP1roG.hpp"


/**
 *  Check the solution of the AP1roG PSEs (with zero initial guess) with reference data from Ayers' implementation
 *  The test system is H2 with HF/6-31G** orbitals
 */
BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // Input the reference data
    const double ref_ap1rog_energy = -1.8696828608304892;
    GQCP::VectorX<double> ref_ap1rog_coefficients (9);
    ref_ap1rog_coefficients << -0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996;


    // Prepare molecular Hamiltonian in the RHF basis
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const auto N_P = h2.numberOfElectrons()/2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Solve the AP1roG PSEs with the initial guess being 0
    const GQCP::AP1roGPSEs pses (sq_hamiltonian, N_P);
    const GQCP::AP1roGPSESolver ap1rog_pse_solver (pses);
    const auto G = ap1rog_pse_solver.solve();  // zero initial guess

    const double electronic_energy = GQCP::QCModel::AP1roG::calculateEnergy(G, sq_hamiltonian);
    const GQCP::VectorX<double> ap1rog_coefficients = G.asVector();


    // Check the results
    BOOST_CHECK(std::abs(electronic_energy - ref_ap1rog_energy) < 1.0e-05);

    for (size_t i = 0; i < 9; i++) {
        BOOST_CHECK(std::abs(ap1rog_coefficients(i) - ref_ap1rog_coefficients(i)) < 1.0e-05);
    }
}


/**
 *  Check the solution of the AP1roG PSEs (with the weak interaction limit as initial guess) with reference data from Ayers' implementation
 *  The test system is H2 with HF/6-31G** orbitals
 */
BOOST_AUTO_TEST_CASE ( h2_631gdp_weak_interaction_limit ) {

    // Input the reference data
    const double ref_ap1rog_energy = -1.8696828608304892;
    GQCP::VectorX<double> ref_ap1rog_coefficients (9);
    ref_ap1rog_coefficients << -0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996;


    // Prepare molecular Hamiltonian in the RHF basis
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2_olsens.xyz");
    const auto N_P = h2.numberOfElectrons() / 2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective (sq_hamiltonian);
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf_parameters.coefficientMatrix());


    // Solve the AP1roG PSEs, with the initial guess being the weak interaction limit coefficients
    const GQCP::AP1roGPSEs pses (sq_hamiltonian, N_P);
    const GQCP::AP1roGPSESolver ap1rog_pse_solver (pses);
    auto G = GQCP::AP1roGGeminalCoefficients::WeakInteractionLimit(sq_hamiltonian, N_P);
    ap1rog_pse_solver.solve(G);  // weak interaction limit coefficients are the initial guess

    const double electronic_energy = GQCP::QCModel::AP1roG::calculateEnergy(G, sq_hamiltonian);
    const GQCP::VectorX<double> ap1rog_coefficients = G.asVector();

    // Check the result
    BOOST_CHECK(std::abs(electronic_energy - ref_ap1rog_energy) < 1.0e-05);

    for (size_t i = 0; i < 9; i++) {
        BOOST_CHECK(std::abs(ap1rog_coefficients(i) - ref_ap1rog_coefficients(i)) < 1.0e-05);
    }
}
