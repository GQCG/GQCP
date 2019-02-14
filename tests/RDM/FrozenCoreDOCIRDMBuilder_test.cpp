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
#define BOOST_TEST_MODULE "FrozenCoreDOCIRDM_test"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>

#include "RDM/RDMCalculator.hpp"
#include "RDM/SelectedRDMBuilder.hpp"

#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FrozenCoreDOCIRDMBuilder.hpp"

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FrozenCoreDOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "utilities/linalg.hpp"


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_one_rdms ) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto ham_par = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    GQCP::FrozenFockSpace fock_space (K, 3, 2);  // dim = 16
    GQCP::SelectedFockSpace selected_fock_space (fock_space);
    GQCP::FrozenCoreDOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the frozen core DOCI and SelectedCI 1-RDMS
    GQCP::SelectedRDMBuilder sci_rdm(selected_fock_space);
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm(fock_space);
    GQCP::OneRDMs one_rdms_s = sci_rdm.calculate1RDMs(coef);
    GQCP::OneRDMs one_rdms = doci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(one_rdms_s.one_rdm.get_matrix_representation().isApprox(one_rdms.one_rdm.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_aa.get_matrix_representation().isApprox(one_rdms.one_rdm_aa.get_matrix_representation()));
    BOOST_CHECK(one_rdms_s.one_rdm_bb.get_matrix_representation().isApprox(one_rdms.one_rdm_bb.get_matrix_representation()));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_two_rdms ) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    auto ham_par = GQCP::HamiltonianParameters::Molecular(H5, "STO-3G");

    GQCP::FrozenFockSpace fock_space (K, 3, 2);  // dim = 16
    GQCP::SelectedFockSpace selected_fock_space (fock_space);
    GQCP::FrozenCoreDOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, ham_par);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    Eigen::VectorXd coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the frozen core DOCI and SelectedCI 1-RDMS
    GQCP::SelectedRDMBuilder sci_rdm(selected_fock_space);
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm(fock_space);
    GQCP::TwoRDMs two_rdms_s = sci_rdm.calculate2RDMs(coef);
    GQCP::TwoRDMs two_rdms = doci_rdm.calculate2RDMs(coef);

    BOOST_CHECK(GQCP::areEqual(two_rdms_s.two_rdm_aaaa.get_matrix_representation(), two_rdms.two_rdm_aaaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(GQCP::areEqual(two_rdms_s.two_rdm_aabb.get_matrix_representation(), two_rdms.two_rdm_aabb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(GQCP::areEqual(two_rdms_s.two_rdm_bbaa.get_matrix_representation(), two_rdms.two_rdm_bbaa.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(GQCP::areEqual(two_rdms_s.two_rdm_bbbb.get_matrix_representation(), two_rdms.two_rdm_bbbb.get_matrix_representation(), 1.0e-06));
    BOOST_CHECK(GQCP::areEqual(two_rdms_s.two_rdm.get_matrix_representation(), two_rdms.two_rdm.get_matrix_representation(), 1.0e-06));
}