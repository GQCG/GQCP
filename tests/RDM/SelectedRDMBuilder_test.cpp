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
#define BOOST_TEST_MODULE "Selected_RDM_test"



#include "RDM/RDMCalculator.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <cpputil.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>


BOOST_AUTO_TEST_CASE ( one_rdms_fci_H2_6_31G ) {

    // Do an H2@FCI//6-31G calculation
    // test if 1-RDM SelectedRDM and FCIRDM are equal
    size_t N_a = 1;
    size_t N_b = 1;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    size_t K = ham_par.get_K();  // 4

    GQCP::ProductFockSpace fock_space (K, N_a, N_b);  // dim = 16
    GQCP::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCP::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the 1-RDM from FCI
    GQCP::RDMCalculator fci_rdm(fock_space);
    GQCP::OneRDMs one_rdms = fci_rdm.calculate1RDMs(coef);


    GQCP::SelectedFockSpace selected_fock_space (fock_space);

    // Get the 1-RDM from SelectedCI
    GQCP::RDMCalculator selected_rdm (selected_fock_space);
    GQCP::OneRDMs one_rdms_s = selected_rdm.calculate1RDMs(coef);


    BOOST_CHECK(one_rdms_s.one_rdm.get_matrix_representation().isApprox(
            one_rdms.one_rdm.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_aa.get_matrix_representation().isApprox(
            one_rdms.one_rdm_aa.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_bb.get_matrix_representation().isApprox(
            one_rdms.one_rdm_bb.get_matrix_representation()));
}



BOOST_AUTO_TEST_CASE ( two_rdms_fci_H2_6_31G ) {

    // Do an H2@FCI//6-31G calculation
    // test if 2-RDM SelectedRDM and FCIRDM are equal
    size_t N_a = 1;
    size_t N_b = 1;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    size_t K = ham_par.get_K();  // 4

    GQCP::ProductFockSpace fock_space (K, N_a, N_b);  // dim = 16
    GQCP::FCI fci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense FCI eigenvalue problem
    GQCP::CISolver ci_solver (fci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the 1-RDM from FCI
    GQCP::RDMCalculator fci_rdm(fock_space);
    GQCP::TwoRDMs two_rdms = fci_rdm.calculate2RDMs(coef);


    GQCP::SelectedFockSpace selected_fock_space (fock_space);

    // Get the 1-RDM from SelectedCI
    GQCP::RDMCalculator selected_rdm (selected_fock_space);
    GQCP::TwoRDMs two_rdms_s = selected_rdm.calculate2RDMs(coef);


    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_aaaa.get_matrix_representation(), two_rdms.two_rdm_aaaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_aabb.get_matrix_representation(), two_rdms.two_rdm_aabb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_bbaa.get_matrix_representation(), two_rdms.two_rdm_bbaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_bbbb.get_matrix_representation(), two_rdms.two_rdm_bbbb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm.get_matrix_representation(), two_rdms.two_rdm.get_matrix_representation(), 1.0e-06));
}


BOOST_AUTO_TEST_CASE ( one_rdms_doci_H2_6_31G ) {

    // Do an H2@doci//6-31G calculation
    // test if 1-RDM SelectedRDM and dociRDM are equal
    size_t N = 1;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    size_t K = ham_par.get_K();  // 4

    GQCP::FockSpace fock_space (K, N);  // dim = 4
    GQCP::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense doci eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the 1-RDM from doci
    GQCP::RDMCalculator doci_rdm(fock_space);
    GQCP::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);


    GQCP::SelectedFockSpace selected_fock_space (fock_space);

    // Get the 1-RDM from SelectedCI
    GQCP::RDMCalculator selected_rdm (selected_fock_space);
    GQCP::OneRDMs one_rdms_s = selected_rdm.calculate1RDMs(coef);


    BOOST_CHECK(one_rdms_s.one_rdm.get_matrix_representation().isApprox(
            one_rdms.one_rdm.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_aa.get_matrix_representation().isApprox(
            one_rdms.one_rdm_aa.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_bb.get_matrix_representation().isApprox(
            one_rdms.one_rdm_bb.get_matrix_representation()));
}



BOOST_AUTO_TEST_CASE ( two_rdms_doci_H2_6_31G ) {

    // Do an H2@doci//6-31G calculation
    // test if 2-RDM SelectedRDM and dociRDM are equal
    size_t N = 1;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    size_t K = ham_par.get_K();  // 4

    GQCP::FockSpace fock_space (K, N);  // dim = 4
    GQCP::DOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense doci eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the 1-RDM from doci
    GQCP::RDMCalculator doci_rdm(fock_space);
    GQCP::TwoRDMs two_rdms = doci_rdm.calculate2RDMs(coef);


    GQCP::SelectedFockSpace selected_fock_space (fock_space);

    // Get the 1-RDM from SelectedCI
    GQCP::RDMCalculator selected_rdm (selected_fock_space);
    GQCP::TwoRDMs two_rdms_s = selected_rdm.calculate2RDMs(coef);


    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_aaaa.get_matrix_representation(), two_rdms.two_rdm_aaaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_aabb.get_matrix_representation(), two_rdms.two_rdm_aabb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_bbaa.get_matrix_representation(), two_rdms.two_rdm_bbaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm_bbbb.get_matrix_representation(), two_rdms.two_rdm_bbbb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(cpputil::linalg::areEqual(two_rdms_s.two_rdm.get_matrix_representation(), two_rdms.two_rdm.get_matrix_representation(), 1.0e-06));
}
