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

#define BOOST_TEST_MODULE "Hubbard"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"


/**
 *  Check if a dense diagonalization using the specialized Hubbard routines produces the same results as when using the unspecialized routines.
 */
BOOST_AUTO_TEST_CASE(Hubbard_specialized_vs_unspecialized_dense_diagonalization) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 4;    // The number of lattice sites.
    const auto N_P = 2;  // The number of electron pairs.

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};

    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Create one dense solver and two environments.
    auto solver = GQCP::EigenproblemSolver::Dense();

    auto specialized_environment = GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);
    auto unspecialized_environment = GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);


    // Optimize the CI model for both environments and check if the energies are equal.
    const auto specialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, specialized_environment).groundStateEnergy();
    const auto unspecialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, unspecialized_environment).groundStateEnergy();
    BOOST_CHECK(std::abs(specialized_energy - unspecialized_energy) < 1.0e-06);
}


/**
 *  Check if a dense diagonalization using the specialized Hubbard routines produces the same results as when using the unspecialized routines, for a larger number of lattice sites.
 */
BOOST_AUTO_TEST_CASE(Hubbard_specialized_vs_unspecialized_dense_diagonalization_large) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 6;    // The number of lattice sites.
    const auto N_P = 3;  // The number of electron pairs.

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};

    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Create one dense solver and two environments.
    auto solver = GQCP::EigenproblemSolver::Dense();

    auto specialized_environment = GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);
    auto unspecialized_environment = GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);


    // Optimize the CI model for both environments and check if the energies are equal.
    const auto specialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, specialized_environment).groundStateEnergy();
    const auto unspecialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, unspecialized_environment).groundStateEnergy();
    BOOST_CHECK(std::abs(specialized_energy - unspecialized_energy) < 1.0e-06);
}


/**
 *  Check if we can reproduce the ground-state energies for a four-site chain from a reference implementation by Ward Poelmans. (https://github.com/wpoely86/Hubbard-GPU)
 */
BOOST_AUTO_TEST_CASE(four_site_chain) {

    // Create the adjacency matrix for a four-site chain.
    const auto K = 4;    // The number of lattice sites.
    const auto N_P = 2;  // The number of electron pairs.
    const auto A = GQCP::AdjacencyMatrix::Linear(K);

    // Set the reference results.
    const double t = 1.0;
    const std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};
    const std::vector<double> E_list {-4.472135955, -3.57536562, -3.202271824, -2.135871608, -1.338959715, -1.004379799, -0.9114974686};


    // Create the appropriate spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    for (size_t i = 0; i < 7; i++) {

        // Create the Hubbard model Hamiltonian.
        const GQCP::HoppingMatrix<double> H {A, t, U_list[i]};
        const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


        // Optimize the CI model, using a dense solver.
        auto solver = GQCP::EigenproblemSolver::Dense();
        auto environment = GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);

        const auto energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();
        BOOST_CHECK(std::abs(energy - E_list[i]) < 1.0e-08);
    }
}


/**
 *  Check if we can reproduce the ground-state energies for a six-site ring from a reference implementation by Ward Poelmans. (https://github.com/wpoely86/Hubbard-GPU)
 */
BOOST_AUTO_TEST_CASE(six_site_ring) {

    // Create the adjacency matrix for a six-site ring.
    const auto K = 6;    // The number of lattice sites.
    const auto N_P = 3;  // The number of electron pairs.
    const auto A = GQCP::AdjacencyMatrix::Cyclic(K);


    // Set the reference results.
    const double t = 1.0;
    const std::vector<double> U_list {0.0, 1.0, 1.5, 3.5, 6.5, 9, 10};
    const std::vector<double> E_list {-8, -6.601158293, -5.978815789, -4.025796251, -2.469458295, -1.836926909, -1.664362733};


    // Create the appropriate spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    for (size_t i = 0; i < 7; i++) {

        // Create the Hubbard model Hamiltonian.
        const GQCP::HoppingMatrix<double> H {A, t, U_list[i]};
        const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};


        // Optimize the CI model, using a dense solver.
        auto solver = GQCP::EigenproblemSolver::Dense();
        auto environment = GQCP::CIEnvironment::Dense(hubbard_hamiltonian, onv_basis);

        const auto energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateEnergy();
        BOOST_CHECK(std::abs(energy - E_list[i]) < 1.0e-08);
    }
}


/**
 *  Check if a Davidson diagonalization using the specialized Hubbard routines produces the same results as when using the unspecialized routines.
 */
BOOST_AUTO_TEST_CASE(Hubbard_specialized_vs_unspecialized_Davidson_diagonalization) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 4;    // The number of lattice sites.
    const auto N_P = 2;  // The number of electron pairs.

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};

    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Create one dense solver and two environments.
    auto solver = GQCP::EigenproblemSolver::Davidson();

    const auto initial_guess = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis).coefficients();
    auto specialized_environment = GQCP::CIEnvironment::Iterative(hubbard_hamiltonian, onv_basis, initial_guess);
    auto unspecialized_environment = GQCP::CIEnvironment::Iterative(hamiltonian, onv_basis, initial_guess);


    // Optimize the CI model for both environments and check if the energies are equal.
    const auto specialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, specialized_environment).groundStateEnergy();
    const auto unspecialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, unspecialized_environment).groundStateEnergy();
    BOOST_CHECK(std::abs(specialized_energy - unspecialized_energy) < 1.0e-06);
}


/**
 *  Check if a Davidson diagonalization using the specialized Hubbard routines produces the same results as when using the unspecialized routines, for a larger number of lattice sites.
 */
BOOST_AUTO_TEST_CASE(Hubbard_specialized_vs_unspecialized_Davidson_diagonalization_large) {

    // Create the Hubbard model Hamiltonian and an appropriate ONV basis.
    const auto K = 6;    // The number of lattice sites.
    const auto N_P = 3;  // The number of electron pairs.

    const auto H = GQCP::HoppingMatrix<double>::Random(K);
    const GQCP::HubbardHamiltonian<double> hubbard_hamiltonian {H};

    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};


    // Create an identical, but 'unspecialized' second-quantized Hamiltonian.
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::FromHubbard(hubbard_hamiltonian);


    // Create one dense solver and two environments.
    auto solver = GQCP::EigenproblemSolver::Davidson();

    const auto initial_guess = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis).coefficients();
    auto specialized_environment = GQCP::CIEnvironment::Iterative(hubbard_hamiltonian, onv_basis, initial_guess);
    auto unspecialized_environment = GQCP::CIEnvironment::Iterative(hamiltonian, onv_basis, initial_guess);


    // Optimize the CI model for both environments and check if the energies are equal.
    const auto specialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, specialized_environment).groundStateEnergy();
    const auto unspecialized_energy = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, unspecialized_environment).groundStateEnergy();
    BOOST_CHECK(std::abs(specialized_energy - unspecialized_energy) < 1.0e-06);
}
