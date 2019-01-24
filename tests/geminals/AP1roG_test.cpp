// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2018  the GQCG developers
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
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "geminals/AP1roG.hpp"

#include "geminals/AP1roGBivariationalSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "properties/expectation_values.hpp"




BOOST_AUTO_TEST_CASE ( calculateOverlap ) {

    size_t N_P = 5;
    size_t K = 8;
    auto N_P_ = static_cast<double>(N_P);
    auto K_ = static_cast<double>(K);


    // Set all bivariational coefficients to 1
    size_t dim = GQCP::AP1roGVariables::numberOfVariables(N_P, K);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(dim);
    for (size_t mu = 0; mu < dim; mu++) {
        q(mu) = 1.0;
    }

    GQCP::BivariationalCoefficients Q {1.0, GQCP::AP1roGVariables(q, N_P, K)};


    // Set the geminal coefficients to i*a
    Eigen::VectorXd g = Eigen::VectorXd::Zero(dim);
    for (size_t i = 0; i < N_P; i++) {
        for (size_t a = N_P; a < K; a++) {
            size_t vector_index = Q.q.vectorIndex(i, a);

            auto i_ = static_cast<double>(i) + 1;
            auto a_ = static_cast<double>(a) + 1;

            g(vector_index) = i_ * a_;
        }
    }

    GQCP::AP1roGGeminalCoefficients G (g, N_P, K);


    // Check with manual formula
    BOOST_CHECK(std::abs(GQCP::calculateOverlap(G, Q) - (1 + 0.25 * N_P_ * (N_P_ + 1) * (K_ * (K_ + 1) - N_P_ * (N_P_ + 1)))) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( energy_as_contraction ) {

    auto h2 = GQCP::Molecule::Readxyz("data/h2_olsens.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G**");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Solve the AP1roG bivariational equations with the initial guess of the geminal coefficients being 0
    GQCP::AP1roGBivariationalSolver bivar_solver (h2, mol_ham_par);
    bivar_solver.solve();
    double electronic_energy = bivar_solver.get_electronic_energy();
    auto G = bivar_solver.get_geminal_coefficients();
    auto Q = bivar_solver.get_bivariational_coefficients();


    // Calculate the 1- and 2-RDM for bivariational AP1roG
    auto D = GQCP::calculate1RDM(G, Q);
    auto d = GQCP::calculate2RDM(G, Q);

    double electronic_energy_by_contraction = GQCP::calculateExpectationValue(mol_ham_par, D, d) - mol_ham_par.get_scalar();  // only the electronic energy
    std::cout << electronic_energy << std::endl;
    std::cout << electronic_energy_by_contraction << std::endl;

    BOOST_CHECK(std::abs(electronic_energy_by_contraction - electronic_energy) < 1.0e-12);
}
