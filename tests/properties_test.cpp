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
#include "LibintCommunicator.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


#include "units.hpp"


BOOST_AUTO_TEST_CASE ( dipole_CO_STO_3G ) {

    // Initialize the molecule and molecular Hamiltonian parameters for CO
    GQCP::Atom C (6, 0.0, 0.0, 0.0);
    GQCP::Atom O (8, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.145));  // from CCCBDB
    std::vector<GQCP::Atom> atoms {C, O};
    GQCP::Molecule CO (atoms);

    auto ao_basis = std::make_shared<GQCP::AOBasis>(CO, "STO-3G");
    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    size_t K = ao_basis->get_number_of_basis_functions();
    size_t N = CO.get_N();

    // Solve the SCF equations
    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, CO);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    double total_energy = rhf.get_electronic_energy() + CO.calculateInternuclearRepulsionEnergy();
    BOOST_REQUIRE(std::abs(total_energy - (-111.225)) < 1.0e-02);  // from CCCBDB, require a correct RHF solution to be found


    // Calculate the RHF 1-RDM in MO basis
    auto D = GQCP::calculateRHF1RDM(K, N);
    auto D_AO = GQCP::calculateRHFAO1RDM(rhf.get_C(), N);

    // Calculate the dipole integrals, and transform them to the MO basis
    auto dipole_components = GQCP::LibintCommunicator::get().calculateDipoleIntegrals(*ao_basis);
    for (auto& dipole_component : dipole_components) {
        dipole_component.transform(rhf.get_C());
    }

    Eigen::Vector3d total_dipole_moment = CO.calculateNuclearDipoleMoment() + GQCP::calculateElectronicDipoleMoment(dipole_components, D);
    BOOST_CHECK(std::abs(total_dipole_moment.norm() - (0.049)) < 1.0e-03);
}
