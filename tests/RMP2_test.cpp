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
#define BOOST_TEST_MODULE "RMP2"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "RMP2.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

BOOST_AUTO_TEST_CASE ( crawdad_sto3g_water ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt)
    double ref_energy_correction = -0.049149636120;


    // Create molecular Hamiltonian parameters in the RHF basis
    auto water = GQCP::Molecule::Readxyz("../tests/data/h2o_crawdad.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(water, "STO-3G");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, water);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(mol_ham_par, water, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( crawdad_sto3g_methane ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/ch4_sto3g/output.txt)
    double ref_energy_correction = -0.056046676165;

    // Create molecular Hamiltonian parameters in the RHF basis
    auto methane = GQCP::Molecule::Readxyz("../tests/data/ch4_crawdad.xyz");
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(methane, "STO-3G");

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, methane);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(mol_ham_par, methane, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}
