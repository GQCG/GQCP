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
#define BOOST_TEST_MODULE "SelectedFockSpace"

#include <boost/test/unit_test.hpp>

#include "FockSpace/SelectedFockSpace.hpp"
#include "WaveFunction/WaveFunctionReader.hpp"


BOOST_AUTO_TEST_CASE ( constructor ) {

    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace (10, 5, 5));

    GQCP::ProductFockSpace fock_space_product (10, 5, 5);
    GQCP::FockSpace fock_space (10, 5);

    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace fock (fock_space_product));
    BOOST_CHECK_NO_THROW(GQCP::SelectedFockSpace fock (fock_space));
}


BOOST_AUTO_TEST_CASE ( addConfiguration ) {

    // Create a faulty expansion: one of the orbitals is different
    GQCP::SelectedFockSpace fock_space (3, 1, 1);

    std::vector<std::string> alpha_set = {"001", "010"};
    std::vector<std::string> beta_set = {"001", "010"};

    BOOST_CHECK_NO_THROW(fock_space.addConfiguration(alpha_set, beta_set));

    // Test throw with one of the sets is not the same size
    std::vector<std::string> beta_set_long = {"001", "010", "100"};
    BOOST_CHECK_THROW(fock_space.addConfiguration(alpha_set, beta_set_long), std::invalid_argument);

    // Test throw with incompatible orbital numbers
    BOOST_CHECK_THROW(fock_space.addConfiguration("0001", "0100"), std::invalid_argument);

    // Test throw with incompatible electron numbers
    BOOST_CHECK_THROW(fock_space.addConfiguration("011", "011"), std::invalid_argument);

    fock_space.addConfiguration(alpha_set, beta_set);


    // Check if the expansions are equal
    // Generate the expected results
    std::string alpha1_ref = "001";
    std::string alpha2_ref = "010";
    std::string beta1_ref = "001";
    std::string beta2_ref = "010";

    // Retrieve the added results
    GQCP::Configuration configuration1 = fock_space.get_configuration(0);
    GQCP::Configuration configuration2 = fock_space.get_configuration(1);

    // Retrieve the string representation of the ONVs
    std::string alpha1_test = configuration1.onv_alpha.asString();
    std::string alpha2_test = configuration2.onv_alpha.asString();
    std::string beta1_test = configuration1.onv_beta.asString();
    std::string beta2_test = configuration2.onv_beta.asString();

    BOOST_CHECK(alpha1_test == alpha1_ref);
    BOOST_CHECK(alpha2_test == alpha2_ref);
    BOOST_CHECK(beta1_test == beta1_ref);
    BOOST_CHECK(beta2_test == beta2_ref);
}


BOOST_AUTO_TEST_CASE ( reader_test ) {

    // We will test if we can construct a selected fock space and a corresponding coefficients
    // through an input file
    GQCP::WaveFunctionReader test_reader ("data/test_GAMESS_expansion");


    // Check read vector is correct
    // Gerenate the expected result
    Eigen::Vector2d test_vector;
    test_vector << 1.0000, 0.0000;

    BOOST_CHECK(test_vector.isApprox(test_reader.get_coefficients()));

    // Check if the expansions are equal
    // Generate the expected results
    std::string alpha1_ref = "0000000000000000000000000000000000000000000001";
    std::string alpha2_ref = "0000000000000000000000000000000000000000000001";
    std::string beta1_ref = "0000000000000000000000000000000000000000000001";
    std::string beta2_ref = "0000000000000000000000000000000000000000000010";

    // Retrieve the read results
    GQCP::Configuration configuration1 = test_reader.get_fock_space().get_configuration(0);
    GQCP::Configuration configuration2 = test_reader.get_fock_space().get_configuration(1);

    // Retrieve the string representation of the ONVs
    std::string alpha1_test = configuration1.onv_alpha.asString();
    std::string alpha2_test = configuration2.onv_alpha.asString();
    std::string beta1_test = configuration1.onv_beta.asString();
    std::string beta2_test = configuration2.onv_beta.asString();

    BOOST_CHECK(alpha1_test == alpha1_ref);
    BOOST_CHECK(alpha2_test == alpha2_ref);
    BOOST_CHECK(beta1_test == beta1_ref);
    BOOST_CHECK(beta2_test == beta2_ref);

}


BOOST_AUTO_TEST_CASE ( Selected_Evaluation_H2O ) {

    // Psi4 and GAMESS' FCI energy for H2O
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian parameters in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(h2o, "STO-3G");
    auto K = mol_ham_par.get_K();

    mol_ham_par.LowdinOrthonormalize();

    GQCP::ProductFockSpace fock_space (K, h2o.numberOfElectrons()/2, h2o.numberOfElectrons()/2);  // dim = 441
    GQCP::SelectedFockSpace selected_fock_space (fock_space);

    GQCP::SquareMatrix<double> hamiltonian = selected_fock_space.evaluateOperatorDense(mol_ham_par, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = selected_fock_space.evaluateOperatorDense(mol_ham_par, false);
    GQCP::VectorX<double> hamiltonian_diagonal = selected_fock_space.evaluateOperatorDiagonal(mol_ham_par);

    // Retrieve lowest eigenvalue (fci solution)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver (hamiltonian);

    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
    double test_energy = self_adjoint_eigensolver.eigenvalues()(0) + internuclear_repulsion_energy;


    BOOST_CHECK(std::abs(test_energy - reference_fci_energy) < 1e-6);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));
}
