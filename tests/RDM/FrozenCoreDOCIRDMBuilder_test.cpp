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

#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FrozenCoreDOCI.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FrozenCoreDOCIRDMBuilder.hpp"
#include "RDM/RDMCalculator.hpp"
#include "RDM/SelectedRDMBuilder.hpp"
#include "Utilities/linalg.hpp"


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_one_rdms ) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (H5, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H5);  // in an AO basis

    GQCP::FrozenFockSpace fock_space (K, 3, 2);
    GQCP::SelectedFockSpace selected_fock_space (fock_space);
    GQCP::FrozenCoreDOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, sq_hamiltonian);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::VectorX<double> coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the frozen core DOCI and SelectedCI 1-RDMS
    GQCP::SelectedRDMBuilder sci_rdm (selected_fock_space);
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm (fock_space);
    GQCP::OneRDMs<double> one_rdms_s = sci_rdm.calculate1RDMs(coef);
    GQCP::OneRDMs<double> one_rdms = doci_rdm.calculate1RDMs(coef);

    BOOST_CHECK(one_rdms_s.one_rdm.isApprox(one_rdms.one_rdm));
    BOOST_CHECK(one_rdms_s.one_rdm_aa.isApprox(one_rdms.one_rdm_aa));
    BOOST_CHECK(one_rdms_s.one_rdm_bb.isApprox(one_rdms.one_rdm_bb));
}


BOOST_AUTO_TEST_CASE ( FrozenCoreDOCI_two_rdms ) {

    size_t K = 5;
    GQCP::Molecule H5 = GQCP::Molecule::HChain(K, 1.1);
    GQCP::RSpinorBasis<double, GQCP::GTOShell> sp_basis (H5, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(sp_basis, H5);  // in an AO basis

    GQCP::FrozenFockSpace fock_space (K, 3, 2);
    GQCP::SelectedFockSpace selected_fock_space (fock_space);
    GQCP::FrozenCoreDOCI doci (fock_space);

    // Specify solver options and solve the eigenvalue problem
    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, sq_hamiltonian);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    GQCP::VectorX<double> coef = ci_solver.get_eigenpair().get_eigenvector();

    // Get the frozen core DOCI and SelectedCI 2-RDMS
    GQCP::SelectedRDMBuilder sci_rdm (selected_fock_space);
    GQCP::FrozenCoreDOCIRDMBuilder doci_rdm (fock_space);
    GQCP::TwoRDMs<double> two_rdms_s = sci_rdm.calculate2RDMs(coef);
    GQCP::TwoRDMs<double> two_rdms = doci_rdm.calculate2RDMs(coef);

    BOOST_CHECK(two_rdms_s.two_rdm_aaaa.isApprox(two_rdms.two_rdm_aaaa, 1.0e-06));
    BOOST_CHECK(two_rdms_s.two_rdm_aabb.isApprox(two_rdms.two_rdm_aabb, 1.0e-06));
    BOOST_CHECK(two_rdms_s.two_rdm_bbaa.isApprox(two_rdms.two_rdm_bbaa, 1.0e-06));
    BOOST_CHECK(two_rdms_s.two_rdm_bbbb.isApprox(two_rdms.two_rdm_bbbb, 1.0e-06));
    BOOST_CHECK(two_rdms_s.two_rdm.isApprox(two_rdms.two_rdm, 1.0e-06));
}
