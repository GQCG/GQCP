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
#define BOOST_TEST_MODULE "DavidsonHubbardSolver"


#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_davidson ) {

    // Check if FCI and Hubbard produce the same results for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(10);

    size_t N = 2;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();

    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 36

    // Create the Hubbard and FCI modules
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::CISolver solver1 (hubbard, mol_ham_par);
    GQCP::CISolver solver2 (fci, mol_ham_par);

    // Solve with Davidson
    Eigen::VectorXd initial_guess = fock_space.randomExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_guess);
    solver1.solve(solver_options);
    solver2.solve(solver_options);

    auto fci_energy = solver2.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = solver1.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( test_Hubbard_vs_FCI_davidson_large ) {

    // Check if FCI and Hubbard produce the same results for Hubbard Hamiltonian parameters

    // Create the Hamiltonian parameters for the triagonal of a Hubbard lattice.
    Eigen::VectorXd triagonal_test = Eigen::VectorXd::Random(21);

    size_t N = 3;
    auto mol_ham_par = GQCP::constructHubbardParameters(triagonal_test);
    auto K = mol_ham_par.get_K();

    GQCP::ProductFockSpace fock_space (K, N, N);  // dim = 400

    // Create the Hubbard and FCI modules
    GQCP::Hubbard hubbard (fock_space);
    GQCP::FCI fci (fock_space);

    GQCP::CISolver solver1 (hubbard, mol_ham_par);
    GQCP::CISolver solver2 (fci, mol_ham_par);

    // Solve with Davidson
    Eigen::VectorXd initial_guess = fock_space.randomExpansion();
    numopt::eigenproblem::DavidsonSolverOptions solver_options (initial_guess);
    solver1.solve(solver_options);
    solver2.solve(solver_options);

    auto fci_energy = solver2.get_eigenpair().get_eigenvalue();
    auto hubbard_energy = solver1.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(fci_energy - (hubbard_energy)) < 1.0e-06);
}
