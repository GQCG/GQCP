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
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain




BOOST_AUTO_TEST_CASE ( Hubbard_h2_sto3g_dense_vs_Davidson ) {

    // Check if the dense Hubbard energy is equal to the Davidson (with matvec) Hubbard energy

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpaceProduct fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 2

    // Create the Hubbard module
    GQCP::Hubbard Hubbard (fock_space);
    GQCP::CISolver ci_solver (Hubbard, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(Hubbard_dense_eigenvalue - Hubbard_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( Hubbard_H2_6_31Gxx_dense_vs_Davidson ) {

    // Check if the dense Hubbard energy is equal to the Davidson (with matvec) Hubbard energy

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31G**");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpaceProduct fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 100

    // Create the Hubbard module
    GQCP::Hubbard Hubbard (fock_space);
    GQCP::CISolver ci_solver (Hubbard, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(Hubbard_dense_eigenvalue - Hubbard_davidson_eigenvalue) < 1.0e-08);
}


BOOST_AUTO_TEST_CASE ( Hubbard_H2O_STO_3G_dense_vs_Davidson ) {

    // Check if the dense Hubbard energy is equal to the Davidson (with matvec) Hubbard energy

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2o ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::FockSpaceProduct fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the Hubbard module
    GQCP::Hubbard Hubbard (fock_space);
    GQCP::CISolver ci_solver (Hubbard, mol_ham_par);

    // Solve Davidson
    Eigen::VectorXd initial_g = fock_space.HartreeFockExpansion();
    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options (initial_g);
    ci_solver.solve(davidson_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto Hubbard_dense_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    BOOST_CHECK(std::abs(Hubbard_dense_eigenvalue - Hubbard_davidson_eigenvalue) < 1.0e-08);
}
