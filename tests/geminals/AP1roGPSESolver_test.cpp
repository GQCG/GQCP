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

#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGPSESolver.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"



BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // We have some reference data from olsens: H2 with HF/6-31G** orbitals
    //  AP1roG energy: -1.8696828608304892
    //  AP1roG coefficients: [-0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996]

    double ref_ap1rog_energy = -1.8696828608304892;
    GQCP::VectorX<double> ref_ap1rog_coefficients (9);
    ref_ap1rog_coefficients << -0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996;


    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("data/h2_olsens.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2, "6-31G**");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters<double>(ao_mol_ham_par, rhf.get_C());
  
  
    // Solve the AP1roG pSE equations with the initial guess being 0
    GQCP::AP1roGPSEs pses (mol_ham_par, h2.get_N()/2);
    GQCP::AP1roGPSESolver ap1rog_pse_solver (pses);
    auto G = ap1rog_pse_solver.solve();
    double electronic_energy = GQCP::calculateAP1roGEnergy(G, mol_ham_par);


    BOOST_CHECK(std::abs(electronic_energy - ref_ap1rog_energy) < 1.0e-05);

    GQCP::VectorX<double> ap1rog_coefficients = G.asVector();
    for (size_t i = 0; i < 9; i++) {
        BOOST_CHECK(std::abs(ap1rog_coefficients(i) - ref_ap1rog_coefficients(i)) < 1.0e-05);
    }
}


BOOST_AUTO_TEST_CASE ( h2_631gdp_weak_interaction_limit ) {

    // We have some reference data from olsens: H2 with HF/6-31G** orbitals
    //  AP1roG energy: -1.8696828608304892
    //  AP1roG coefficients: [-0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996]

    double ref_ap1rog_energy = -1.8696828608304892;
    GQCP::VectorX<double> ref_ap1rog_coefficients (9);
    ref_ap1rog_coefficients << -0.05949796, -0.05454253, -0.03709503, -0.02899231, -0.02899231, -0.01317386, -0.00852702, -0.00852702, -0.00517996;


    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("data/h2_olsens.xyz");
    size_t N_P = h2.get_N() / 2;
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2, "6-31G**");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters<double>(ao_mol_ham_par, rhf.get_C());


    // Solve the AP1roG pSE equations, with the initial guess being the weak interaction limit coefficients
    GQCP::AP1roGPSEs pses (mol_ham_par, h2.get_N()/2);
    GQCP::AP1roGPSESolver ap1rog_pse_solver (pses);
    auto G_initial = GQCP::AP1roGGeminalCoefficients::WeakInteractionLimit(mol_ham_par, N_P);
    auto G = ap1rog_pse_solver.solve();
    double electronic_energy = GQCP::calculateAP1roGEnergy(G, mol_ham_par);


    BOOST_CHECK(std::abs(electronic_energy - ref_ap1rog_energy) < 1.0e-05);

    GQCP::VectorX<double> ap1rog_coefficients = G.asVector();
    for (size_t i = 0; i < 9; i++) {
        BOOST_CHECK(std::abs(ap1rog_coefficients(i) - ref_ap1rog_coefficients(i)) < 1.0e-05);
    }
}
