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

#define BOOST_TEST_MODULE "FCI_RDM_test"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedDMCalculator.hpp"
#include "Processing/DensityMatrices/CIDMCalculators/SpinResolvedSelectedDMCalculator.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/**
 *  Check if the trace of the spin-resolved 1-DMs yields the appropriate number of electrons.
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441.
 */
BOOST_AUTO_TEST_CASE(density_matrices_traces) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = 5;
    const size_t N_beta = 5;
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs, calculate the traces and check if they match the expected result.
    const GQCP::SpinResolvedDMCalculator spin_resolved_rdm_builder {onv_basis};
    const auto one_rdms = spin_resolved_rdm_builder.calculate1RDMs(linear_expansion.coefficients());

    BOOST_CHECK(std::abs(one_rdms.alpha().trace() - N_alpha) < 1.0e-12);
    BOOST_CHECK(std::abs(one_rdms.beta().trace() - N_beta) < 1.0e-12);
    BOOST_CHECK(std::abs(one_rdms.spinSummed().trace() - (N_alpha + N_beta)) < 1.0e-12);


    // Calculate the 2-DMs, calculate the traces and check if they match the expected result.
    const auto two_rdms = spin_resolved_rdm_builder.calculate2RDMs(linear_expansion.coefficients());

    BOOST_CHECK(std::abs(two_rdms.two_rdm_aaaa.trace() - N_alpha * (N_alpha - 1)) < 1.0e-12);
    BOOST_CHECK(std::abs(two_rdms.two_rdm_aabb.trace() - N_alpha * N_beta) < 1.0e-12);
    BOOST_CHECK(std::abs(two_rdms.two_rdm_bbaa.trace() - N_beta * N_alpha) < 1.0e-12);
    BOOST_CHECK(std::abs(two_rdms.two_rdm_bbbb.trace() - N_beta * (N_beta - 1)) < 1.0e-12);


    GQCP::OneDM<double> D_from_reduction = (1.0 / (N - 1)) * two_rdms.two_rdm.reduce();
    BOOST_CHECK(one_rdms.spinSummed().isApprox(D_from_reduction, 1.0e-12));
}


/**
 *  Check if the 1-DM can be calculated from a reduction of the 2-DM.
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441.
 */
BOOST_AUTO_TEST_CASE(One_DM_from_2_DM) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = 5;
    const size_t N_beta = 5;
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1 and 2-DMs, calculate the traces and check if they match the expected result.
    const GQCP::SpinResolvedDMCalculator spin_resolved_rdm_builder {onv_basis};
    const auto D = spin_resolved_rdm_builder.calculate1RDMs(linear_expansion.coefficients()).spinSummed();
    const auto d = spin_resolved_rdm_builder.calculate2RDMs(linear_expansion.coefficients()).two_rdm;


    const auto D_from_reduction = (1.0 / (N - 1)) * d.reduce();
    BOOST_CHECK(D_from_reduction.isApprox(D, 1.0e-12));
}


/**
 *  Check if the energy is equal to the expectation value of the Hamiltonian.
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441.
 */
BOOST_AUTO_TEST_CASE(energy_expectation_value_Hamiltonian) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = 5;
    const size_t N_beta = 5;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto qc_structure = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment);
    const auto energy = qc_structure.groundStateEnergy();
    const auto& linear_expansion = qc_structure.groundStateParameters();


    // Calculate the 1 and 2-DMs and calculate the expectation value of the Hamiltonian.
    const GQCP::SpinResolvedDMCalculator spin_resolved_rdm_builder {onv_basis};
    const auto D = spin_resolved_rdm_builder.calculate1RDMs(linear_expansion.coefficients()).spinSummed();
    const auto d = spin_resolved_rdm_builder.calculate2RDMs(linear_expansion.coefficients()).two_rdm;

    const auto energy_by_contraction = sq_hamiltonian.calculateExpectationValue(D, d);
    BOOST_CHECK(std::abs(energy - energy_by_contraction) < 1.0e-12);
}


/**
 *  Check if the 1- and 2-DMs for a full spin-resolved ONV basis are equal to the 'selected' case.
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441. However, we're choosing a different number of alpha and beta electrons. (N_alpha = 4, N_beta = 6)
 */
BOOST_AUTO_TEST_CASE(specialized_vs_selected_DMs) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = 4;
    const size_t N_beta = 6;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const GQCP::SpinResolvedDMCalculator spin_resolved_rdm_builder {onv_basis};
    const auto one_rdms_specialized = spin_resolved_rdm_builder.calculate1RDMs(linear_expansion.coefficients());

    const GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};
    const GQCP::SpinResolvedSelectedDMCalculator selected_rdm_builder {selected_onv_basis};
    const auto one_rdms_selected = selected_rdm_builder.calculate1RDMs(linear_expansion.coefficients());

    BOOST_CHECK(one_rdms_specialized.spinSummed().isApprox(one_rdms_selected.spinSummed(), 1.0e-12));
    BOOST_CHECK(one_rdms_specialized.alpha().isApprox(one_rdms_selected.alpha(), 1.0e-12));
    BOOST_CHECK(one_rdms_specialized.beta().isApprox(one_rdms_selected.beta(), 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto two_rdms_specialized = spin_resolved_rdm_builder.calculate2RDMs(linear_expansion.coefficients());
    const auto two_rdms_selected = selected_rdm_builder.calculate2RDMs(linear_expansion.coefficients());

    BOOST_CHECK(two_rdms_specialized.two_rdm_aaaa.isApprox(two_rdms_selected.two_rdm_aaaa, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_aabb.isApprox(two_rdms_selected.two_rdm_aabb, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_bbaa.isApprox(two_rdms_selected.two_rdm_bbaa, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm_bbbb.isApprox(two_rdms_selected.two_rdm_bbbb, 1.0e-12));
    BOOST_CHECK(two_rdms_specialized.two_rdm.isApprox(two_rdms_selected.two_rdm, 1.0e-12));
}
