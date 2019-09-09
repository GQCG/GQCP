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
#define BOOST_TEST_MODULE "RMP2"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RMP2.hpp"


BOOST_AUTO_TEST_CASE ( crawdad_sto3g_water ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/h2o_sto3g/output.txt)
    double ref_energy_correction = -0.049149636120;


    // Create the molecular Hamiltonian in the RHF basis
    auto water = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (water, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, water);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, water);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    sq_hamiltonian.transform(rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(sq_hamiltonian, water, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( crawdad_sto3g_methane ) {

    // Get the reference data from crawdad (http://sirius.chem.vt.edu/~crawdad/programming/project4/ch4_sto3g/output.txt)
    double ref_energy_correction = -0.056046676165;

    // Create the molecular Hamiltonian in the RHF basis
    auto methane = GQCP::Molecule::ReadXYZ("data/ch4_crawdad.xyz");
    GQCP::SingleParticleBasis<double, GQCP::GTOShell> sp_basis (methane, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, methane);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, methane);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    sq_hamiltonian.transform(rhf.get_C());


    // Check if the RMP2 correction is correct
    double energy_correction = GQCP::calculateRMP2EnergyCorrection(sq_hamiltonian, methane, rhf);
    BOOST_CHECK(std::abs(energy_correction - ref_energy_correction) < 1.0e-08);
}
