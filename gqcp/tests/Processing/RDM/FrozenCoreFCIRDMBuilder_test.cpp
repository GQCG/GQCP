// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE "FrozenCoreFCIRDM_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/FrozenCoreFCIRDMBuilder.hpp"
#include "Processing/RDM/SelectedRDMBuilder.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/**
 *  Check if the 1- and 2-DMs for a frozen core spin-resolved ONV basis are equal to the 'selected' case.
 *  The system of interest is a linear chain of 5 H atoms, 1.1 bohr apart, using an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(FrozenCoreFCI_one_rdms) {

    // Set up the molecular Hamiltonian for H5//STO-3G in the Löwdin basis.
    const GQCP::Molecule molecule = GQCP::Molecule::HChain(5, 1.1);
    const auto N_alpha = 3;
    const auto N_beta = 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense frozen core CI calculation for 2 frozen orbitals.
    const GQCP::SpinResolvedFrozenONVBasis onv_basis {K, N_alpha, N_beta, 2};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedFrozenONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const GQCP::FrozenCoreFCIRDMBuilder spin_resolved_rdm_builder {onv_basis};
    const auto one_rdms_specialized = spin_resolved_rdm_builder.calculate1RDMs(linear_expansion.coefficients());

    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};
    const GQCP::SelectedRDMBuilder selected_rdm_builder {selected_onv_basis};
    const auto one_rdms_selected = selected_rdm_builder.calculate1RDMs(linear_expansion.coefficients());

    BOOST_CHECK(one_rdms_specialized.one_rdm.isApprox(one_rdms_selected.one_rdm, 1.0e-12));
    BOOST_CHECK(one_rdms_specialized.one_rdm_aa.isApprox(one_rdms_selected.one_rdm_aa, 1.0e-12));
    BOOST_CHECK(one_rdms_specialized.one_rdm_bb.isApprox(one_rdms_selected.one_rdm_bb, 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto two_rdms_specialized = spin_resolved_rdm_builder.calculate2RDMs(linear_expansion.coefficients());
    const auto two_rdms_selected = selected_rdm_builder.calculate2RDMs(linear_expansion.coefficients());

    BOOST_CHECK(two_rdms_specialized.two_rdm_aaaa.isApprox(two_rdms_selected.two_rdm_aaaa, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_aabb.isApprox(two_rdms_selected.two_rdm_aabb, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_bbaa.isApprox(two_rdms_selected.two_rdm_bbaa, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_bbbb.isApprox(two_rdms_selected.two_rdm_bbbb, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm.isApprox(two_rdms_selected.two_rdm, 1.0e-12));
}
