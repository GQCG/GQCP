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
#define BOOST_TEST_MODULE "atomic_decomposition"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "RHF/PlainRHFSCFSolver.hpp"
#include "HamiltonianParameters/AtomicDecompositionParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "RDM/RDMCalculator.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "properties/expectation_values.hpp"
#include "units.hpp"



BOOST_AUTO_TEST_CASE ( decomposition_BeH_cation_STO_3G ) {

    // Create the molecular Hamiltonian parameters in an AO basis
    GQCP::Atom Be(4, 0.0, 0.0, 0.0);
    GQCP::Atom H(1, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms{Be, H};
    GQCP::Molecule BeH(atoms, +1);
    GQCP::AtomicDecompositionParameters adp(BeH, "STO-3G");
    auto mol_ham_par = adp.molecular_hamiltonian_parameters;
    auto K = mol_ham_par.get_K();
    double repulsion = BeH.calculateInternuclearRepulsionEnergy();
    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver(mol_ham_par, BeH);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    const auto &T = rhf.get_C();

    // Transform the ham_par
    mol_ham_par.transform(T);

    // Create the FCI module
    GQCP::ProductFockSpace fock_space(K, BeH.get_N() / 2, BeH.get_N() / 2);  // dim = 441
    GQCP::FCI fci(fock_space);
    GQCP::CISolver ci_solver(fci, mol_ham_par);

    // Solve Davidson
    GQCP::DavidsonSolverOptions solver_options(fock_space.HartreeFockExpansion());
    ci_solver.solve(solver_options);

    // Retrieve the eigenpair
    auto fci_coeff = ci_solver.get_eigenpair().get_eigenvector();
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Decomposition calculations require RDM-tracing
    GQCP::RDMCalculator rdm_calc(fock_space);
    rdm_calc.set_coefficients(fci_coeff);

    auto one_rdm = rdm_calc.calculate1RDMs().one_rdm;
    auto two_rdm = rdm_calc.calculate2RDMs().two_rdm;

    // Transform rdms to the AObasis
    GQCP::OneRDM<double> ao_one_rdm = T * (one_rdm) * T.adjoint();
    GQCP::TwoRDM<double> ao_two_rdm = two_rdm;
    ao_two_rdm.fourModeMultiplication<double>(T.adjoint().conjugate(), T.adjoint(), T.adjoint().conjugate(), T.adjoint());

    double self_energy_a = GQCP::calculateExpectationValue(adp.net_atomic_parameters[0], ao_one_rdm, ao_two_rdm);
    double self_energy_b = GQCP::calculateExpectationValue(adp.net_atomic_parameters[1], ao_one_rdm, ao_two_rdm);
    double interaction_energy_ab = GQCP::calculateExpectationValue(adp.interaction_parameters[0], ao_one_rdm,
                                                                   ao_two_rdm);
    double total_energy_a = GQCP::calculateExpectationValue(adp.fragment_parameters[0], ao_one_rdm, ao_two_rdm);
    double total_energy_b = GQCP::calculateExpectationValue(adp.fragment_parameters[1], ao_one_rdm, ao_two_rdm);

    BOOST_CHECK(std::abs(total_energy_a + total_energy_b - fci_energy) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_a + self_energy_b + interaction_energy_ab - fci_energy - repulsion) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_a + interaction_energy_ab / 2 - total_energy_a - repulsion / 2) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_b + interaction_energy_ab / 2 - total_energy_b - repulsion / 2) < 1.0e-010);

}



BOOST_AUTO_TEST_CASE ( decomposition_BeH_cation_STO_3G_Mario ) {

    // Create the molecular Hamiltonian parameters in an AO basis
    GQCP::Atom Be(4, 0.0, 0.0, 0.0);
    GQCP::Atom H(1, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134));  // from CCCBDB, STO-3G geometry
    std::vector<GQCP::Atom> atoms{Be, H};
    GQCP::Molecule BeH(atoms, +1);
    GQCP::AtomicDecompositionParameters adp = GQCP::AtomicDecompositionParameters::Mario(BeH, "STO-3G");
    auto mol_ham_par = adp.molecular_hamiltonian_parameters;
    auto K = mol_ham_par.get_K();
    double repulsion = BeH.calculateInternuclearRepulsionEnergy();
    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver(mol_ham_par, BeH);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    const auto &T = rhf.get_C();

    // Transform the ham_par
    mol_ham_par.transform(T);

    // Create the FCI module
    GQCP::ProductFockSpace fock_space(K, BeH.get_N() / 2, BeH.get_N() / 2);  // dim = 441
    GQCP::FCI fci(fock_space);
    GQCP::CISolver ci_solver(fci, mol_ham_par);

    // Solve Davidson
    GQCP::DavidsonSolverOptions solver_options(fock_space.HartreeFockExpansion());
    ci_solver.solve(solver_options);

    // Retrieve the eigenpair
    auto fci_coeff = ci_solver.get_eigenpair().get_eigenvector();
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Decomposition calculations require RDM-tracing
    GQCP::RDMCalculator rdm_calc(fock_space);
    rdm_calc.set_coefficients(fci_coeff);

    auto one_rdm = rdm_calc.calculate1RDMs().one_rdm;
    auto two_rdm = rdm_calc.calculate2RDMs().two_rdm;

    // Transform rdms to the AObasis
    GQCP::OneRDM<double> ao_one_rdm = T * (one_rdm) * T.adjoint();
    GQCP::TwoRDM<double> ao_two_rdm = two_rdm;
    ao_two_rdm.fourModeMultiplication<double>(T.adjoint().conjugate(), T.adjoint(), T.adjoint().conjugate(), T.adjoint());

    double self_energy_a = GQCP::calculateExpectationValue(adp.net_atomic_parameters[0], ao_one_rdm, ao_two_rdm);
    double self_energy_b = GQCP::calculateExpectationValue(adp.net_atomic_parameters[1], ao_one_rdm, ao_two_rdm);
    double interaction_energy_ab = GQCP::calculateExpectationValue(adp.interaction_parameters[0], ao_one_rdm,
                                                                   ao_two_rdm);
    double total_energy_a = GQCP::calculateExpectationValue(adp.fragment_parameters[0], ao_one_rdm, ao_two_rdm);
    double total_energy_b = GQCP::calculateExpectationValue(adp.fragment_parameters[1], ao_one_rdm, ao_two_rdm);

    BOOST_CHECK(std::abs(total_energy_a + total_energy_b - fci_energy) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_a + self_energy_b + interaction_energy_ab - fci_energy - repulsion) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_a + interaction_energy_ab / 2 - total_energy_a - repulsion / 2) < 1.0e-010);
    BOOST_CHECK(std::abs(self_energy_b + interaction_energy_ab / 2 - total_energy_b - repulsion / 2) < 1.0e-010);

}
