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
#define BOOST_TEST_MODULE "DenseDOCISolver"


#include "CISolver/CISolver.hpp"
#include "FockSpace/FockSpaceProduct.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_Hubbard ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

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

    GQCP::FockSpaceProduct fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 2

    // Create the Hubbard module
    GQCP::Hubbard Hubbard (fock_space);

    Eigen::VectorXd diagonal1 = Hubbard.calculateDiagonal(mol_ham_par);

    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

    // Rotate the hampar using the random unitary matrix
    mol_ham_par.rotate(U_random);

    Eigen::VectorXd diagonal2 = Hubbard.calculateDiagonal(mol_ham_par);

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}

