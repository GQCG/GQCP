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
#define BOOST_TEST_MODULE "properties"

#include "properties.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( dipole ) {

    // Initialize the molecule and molecular Hamiltonian parameters for CO (interatomic distance from Szabo, p201)
    GQCP::Atom C (6, 0.0, 0.0, 0.0);
    GQCP::Atom O (8, 0.0, 0.0, 2.166);
    std::vector<GQCP::Atom> atoms {C, O};
    GQCP::Molecule CO (atoms);

    auto ao_basis = std::make_shared<GQCP::AOBasis>(CO, "STO-3G");
    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    // Solve the SCF equations
    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, CO);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();


    // Calculate the RHF 1-RDM in MO basis

    // Calculate the dipole integrals, and transform them to the MO basis

}
