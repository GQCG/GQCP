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

#define BOOST_TEST_MODULE "FCI"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/*
 *  MARK: Dense DOCI calculations
 */

/**
 *  Check if we can reproduce the DOCI energy for BeH+, using a dense solver. The dimension of the seniority zero sub Fock space is 120.
 *  The reference values are obtained from Caitlin Lanssens.
 */
BOOST_AUTO_TEST_CASE(DOCI_BeH_cation_dense) {

    const double reference_energy = -14.8782216937;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/beh_cation_631g_caitlin.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 4 electrons, so 2 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 2};


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 1.5900757460937498e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for LiH, using a dense solver. The dimension of the seniority zero sub Fock space is 120.
 *  The reference values are obtained from Caitlin Lanssens.
 */
BOOST_AUTO_TEST_CASE(DOCI_LiH_dense) {

    const double reference_energy = -8.0029560313;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/lih_631g_caitlin.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 4 electrons, so 2 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 2};


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 9.6074293445896852e-01;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for Li2, using a dense solver. The dimension of the seniority zero sub Fock space is 816.
 *  The reference values are obtained from Caitlin Lanssens.
 */
BOOST_AUTO_TEST_CASE(DOCI_Li2_dense) {

    const double reference_energy = -15.1153976060;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/li2_321g_klaas.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 6 electrons, so 3 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 3};


    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 3.0036546888874875e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/*
 *  MARK: Davidson DOCI calculations
 */


/**
 *  Check if we can reproduce the DOCI energy for H2O//STO-3G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 21.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_h2o_sto3g_Davidson) {

    const double reference_energy = -74.9671366903;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/h2o_sto3g_klaas.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 10 electrons, so 5 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 5};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 9.7794061444134091e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for BeH+//6-31G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 120.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_BeH_cation_Davidson) {

    const double reference_energy = -14.8782216937;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/beh_cation_631g_caitlin.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 4 electrons, so 2 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 2};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 1.5900757460937498e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for N2//STO-3G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 120.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_N2_STO_3G_Davidson) {

    const double reference_energy = -107.5813316864;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/n2_sto-3g_klaas.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 14 electrons, so 7 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 7};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 2.3786407766990290e+01;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for LiH//6-31G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 120.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_LiH_6_31G_Davidson) {

    const double reference_energy = -8.0029560313;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/lih_631g_caitlin.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 4 electrons, so 2 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 2};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 9.6074293445896852e-01;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for Li2//3-21G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 816.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_Li2_3_21G_Davidson) {

    const double reference_energy = -15.1153976060;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/li2_321g_klaas.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 6 electrons, so 3 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 3};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 3.0036546888874875e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}


/**
 *  Check if we can reproduce the DOCI energy for Li2//3-21G, using a Davidson solver. The dimension of the seniority zero sub Fock space is 1287.
 *  The reference values are obtained from Klaas Gunst.
 */
BOOST_AUTO_TEST_CASE(DOCI_H2O_6_31G_Davidson) {

    const double reference_energy = -76.0125161011;

    // Read in the molecular Hamiltonian from a FCIDUMP file.
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/h2o_631g_klaas.FCIDUMP");
    const auto K = sq_hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    // The species contains 10 electrons, so 5 electron pairs.
    // Construct an appropriate seniority-zero ONV basis.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, 5};


    // Create a Davidson solver and corresponding environment and put them together in the QCMethod.
    const auto initial_guess = GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>::HartreeFock(onv_basis).coefficients();
    auto environment = GQCP::CIEnvironment::Iterative(sq_hamiltonian, onv_basis, initial_guess);
    auto solver = GQCP::EigenproblemSolver::Davidson();
    const auto electronic_energy = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();


    // Check our result with the reference.
    const double internuclear_repulsion_energy = 9.7794061444134091e+00;
    const auto energy = electronic_energy + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(energy - (reference_energy)) < 1.0e-09);
}
