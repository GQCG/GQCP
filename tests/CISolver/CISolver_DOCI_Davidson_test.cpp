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
#define BOOST_TEST_MODULE "DavidsonDOCISolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( DOCI_h2o_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O@STO-3G
    double reference_doci_energy = -74.9671366903;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/h2o_sto3g_klaas.FCIDUMP");

    // The species contains 10 electrons and 7 basis functions, this requires a single Fock Space of 7 orbitals and 5 electrons
    GQCP::FockSpace fock_space (ham_par.get_K(), 5);  // dim = 21

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_h2_sto3g_dense_vs_Davidson ) {

    // Create the molecular Hamiltonian parameters in an AO basis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "STO-3G");

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpace fock_space (mol_ham_par.get_K(), h2.get_N()/2);  // dim = 2

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);
    GQCP::CISolver ci_solver (doci, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto doci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto doci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(doci_dense_eigenvalue - doci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( DOCI_h2_631g_dense_vs_Davidson ) {

    // Create the molecular Hamiltonian parameters in an AO basis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2, "6-31G");

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpace fock_space (mol_ham_par.get_K(), h2.get_N()/2);  // dim = 4

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);
    GQCP::CISolver ci_solver (doci, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto doci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto doci_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(doci_dense_eigenvalue - doci_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for BeH+
    double reference_doci_energy = -14.8782216937;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/beh_cation_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCP::FockSpace fock_space (ham_par.get_K(), 2);  // dim = 120

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_n2_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for N2
    double reference_doci_energy = -107.5813316864;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/n2_sto-3g_klaas.FCIDUMP");

    GQCP::FockSpace fock_space (ham_par.get_K(), 7);  // dim = 120

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  2.3786407766990290E+01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for LiH
    double reference_doci_energy = -8.0029560313;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/lih_631g_caitlin.FCIDUMP");

    GQCP::FockSpace fock_space (ham_par.get_K(), 2);  // dim = 120

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_li2_321g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for Li2
    double reference_doci_energy = -15.1153976060;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/li2_321g_klaas.FCIDUMP");

    GQCP::FockSpace fock_space (ham_par.get_K(), 3);  // dim = 816

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_h2o_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O
    double reference_doci_energy = -76.0125161011;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::HamiltonianParameters::ReadFCIDUMP("../tests/data/h2o_631g_klaas.FCIDUMP");

    GQCP::FockSpace fock_space (ham_par.get_K(), 5);  // dim = 1287

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}

/*
BOOST_AUTO_TEST_CASE ( DOCI_lif_631g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for LiF
    double reference_doci_energy = -107.0007150075;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCP::readFCIDUMPFile("../tests/data/lif_631g_klaas.FCIDUMP");

    GQCP::FockSpace fock_space (ham_par.get_K(), 6);  // dim = 376740

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the Davidson DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_g);
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =   9.1249103487674024e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
*/
