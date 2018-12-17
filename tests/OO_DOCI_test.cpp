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
#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "DOCINewtonOrbitalOptimizer.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RDM/FCIRDMBuilder.hpp"
#include "CISolver/CISolver.hpp"

// dim = 2 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_sto_3g ) {

    // Check if OO-DOCI = FCI for a two-electron system
    double reference_fci_energy = -1.13726333769813;


    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("../tests/data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "STO-3G");
    auto K = ao_mol_ham_par.get_K();

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Do the DOCI orbital optimization using specified solver options
    GQCP::FockSpace fock_space (K, h2.get_N()/2);  // dim = 120
    GQCP::DOCI doci (fock_space);
    GQCP::DenseSolverOptions solver_options;
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);
    orbital_optimizer.solve(solver_options);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 4 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31g ) {

    // Check if OO-DOCI = FCI for a two-electron system, starting from the FCI naturals
    double reference_fci_energy = -1.15168629203274;


    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("../tests/data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G");
    auto K = ao_mol_ham_par.get_K();

    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Transform the molecular Hamiltonian parameters to the FCI natural basis
    size_t N_a = h2.get_N() / 2;
    size_t N_b = h2.get_N() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, mol_ham_par);
    GQCP::DenseSolverOptions solver_options;
    fci_solver.solve(solver_options);

    Eigen::VectorXd coef = fci_solver.get_wavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;
    Eigen::MatrixXd U = one_rdm.diagonalize();

    mol_ham_par.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.get_N()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);
    orbital_optimizer.solve(solver_options);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("../tests/data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G**");
    auto K = ao_mol_ham_par.get_K();


    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Transform the molecular Hamiltonian parameters to the FCI natural basis
    size_t N_a = h2.get_N() / 2;
    size_t N_b = h2.get_N() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, mol_ham_par);
    GQCP::DenseSolverOptions solver_options;
    fci_solver.solve(solver_options);

    Eigen::VectorXd coef = fci_solver.get_wavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;
    Eigen::MatrixXd U = one_rdm.diagonalize();

    mol_ham_par.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.get_N()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);
    orbital_optimizer.solve(solver_options);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx_Davidson ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare molecular Hamiltonian parameters in the RHF basis
    auto h2 = GQCP::Molecule::Readxyz("../tests/data/h2_cristina.xyz");
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();  // 0.713176780299327
    auto ao_mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G**");
    auto K = ao_mol_ham_par.get_K();


    GQCP::PlainRHFSCFSolver plain_scf_solver (ao_mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Transform the molecular Hamiltonian parameters to the FCI natural basis
    size_t N_a = h2.get_N() / 2;
    size_t N_b = h2.get_N() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, mol_ham_par);
    GQCP::DenseSolverOptions solver_options;
    fci_solver.solve(solver_options);

    Eigen::VectorXd coef = fci_solver.get_wavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;
    Eigen::MatrixXd U = one_rdm.diagonalize();

    mol_ham_par.rotate(U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.get_N()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    Eigen::VectorXd initial_g = doci_fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);
    orbital_optimizer.solve(davidson_solver_options);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}
