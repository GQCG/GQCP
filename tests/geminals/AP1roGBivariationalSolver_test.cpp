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
#define BOOST_TEST_MODULE "AP1roGBivariationalSolver"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "geminals/AP1roGBivariationalSolver.hpp"

#include "RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( h2_631gdp ) {

    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("data/h2_olsens.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G**");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Solve the AP1roG bivariational equations with the initial guess being 0 (there are no reference data)
    GQCP::AP1roGBivariationalSolver ap1rog_bivar_solver (h2, mol_ham_par);
    ap1rog_bivar_solver.solve();

    double electronic_energy = ap1rog_bivar_solver.get_electronic_energy();
    Eigen::VectorXd ap1rog_coefficients = ap1rog_bivar_solver.get_geminal_coefficients().asVector();
    auto bivar_coeff = ap1rog_bivar_solver.get_bivariational_coefficients();
}
