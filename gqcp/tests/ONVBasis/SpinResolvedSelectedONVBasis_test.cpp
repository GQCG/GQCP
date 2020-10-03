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

#define BOOST_TEST_MODULE "SpinResolvedSelectedONVBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"


/**
 *  Test the general functionality of the addONV function, by testing throws and retrieving configurations.
 */
BOOST_AUTO_TEST_CASE(addONV) {

    // Create a faulty expansion: one of the orbitals is different
    GQCP::SpinResolvedSelectedONVBasis fock_space {3, 1, 1};

    std::vector<std::string> alpha_set {"001", "010"};
    std::vector<std::string> beta_set {"001", "010"};

    BOOST_CHECK_NO_THROW(fock_space.addONV(alpha_set, beta_set));

    // Test throw with one of the sets is not the same size
    std::vector<std::string> beta_set_long = {"001", "010", "100"};
    BOOST_CHECK_THROW(fock_space.addONV(alpha_set, beta_set_long), std::invalid_argument);

    // Test throw with incompatible orbital numbers
    BOOST_CHECK_THROW(fock_space.addONV("0001", "0100"), std::invalid_argument);

    // Test throw with incompatible electron numbers
    BOOST_CHECK_THROW(fock_space.addONV("011", "011"), std::invalid_argument);

    fock_space.addONV(alpha_set, beta_set);


    // Check if the expansions are equal
    // Generate the expected results
    std::string alpha1_ref = "001";
    std::string alpha2_ref = "010";
    std::string beta1_ref = "001";
    std::string beta2_ref = "010";

    // Retrieve the added results
    GQCP::SpinResolvedONV configuration1 = fock_space.onvWithIndex(0);
    GQCP::SpinResolvedONV configuration2 = fock_space.onvWithIndex(1);

    // Retrieve the string representation of the ONVs
    std::string alpha1_test = configuration1.onv(GQCP::Spin::alpha).asString();
    std::string alpha2_test = configuration2.onv(GQCP::Spin::alpha).asString();
    std::string beta1_test = configuration1.onv(GQCP::Spin::beta).asString();
    std::string beta2_test = configuration2.onv(GQCP::Spin::beta).asString();

    BOOST_CHECK(alpha1_test == alpha1_ref);
    BOOST_CHECK(alpha2_test == alpha2_ref);
    BOOST_CHECK(beta1_test == beta1_ref);
    BOOST_CHECK(beta2_test == beta2_ref);
}


/**
 *  Evaluate the Hamiltonian in a selected ONV basis in which all configurations are selected (Full CI)
 *  Compare the evaluation of a direct matrix vector product to that of the matrix vector product evaluations and test the lowest eigenvalue against of the evaluated Hamiltonian a reference value
 */
BOOST_AUTO_TEST_CASE(Selected_Evaluation_H2O) {

    // Psi4 and GAMESS' FCI energy for H2O
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {h2o, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in the Löwdin basis
    auto K = sq_hamiltonian.numberOfOrbitals();

    GQCP::SpinResolvedONVBasis fock_space {K, h2o.numberOfElectrons() / 2, h2o.numberOfElectrons() / 2};  // dim = 441
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {fock_space};


    // Evaluate the dense Hamiltonian
    GQCP::SquareMatrix<double> hamiltonian = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = selected_fock_space.evaluateOperatorDense(sq_hamiltonian, false);

    // Evaluate the diagonal
    GQCP::VectorX<double> hamiltonian_diagonal = selected_fock_space.evaluateOperatorDiagonal(sq_hamiltonian);

    // Evaluate the matvec
    GQCP::VectorX<double> matvec_evaluation = selected_fock_space.evaluateOperatorMatrixVectorProduct(sq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);

    // Calculate the explicit matvec with the dense evaluations
    GQCP::VectorX<double> matvec_reference = hamiltonian * hamiltonian_diagonal;

    // Retrieve the lowest eigenvalue (FCI solution)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver {hamiltonian};

    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
    double test_energy = self_adjoint_eigensolver.eigenvalues()(0) + internuclear_repulsion_energy;

    // Test the energy with the reference
    BOOST_CHECK(std::abs(test_energy - reference_fci_energy) < 1e-6);

    // Test if the non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));

    // Test if the matvecs are identical
    BOOST_CHECK(matvec_evaluation.isApprox(matvec_reference));
}


/**
 *  Evaluate the an unrestricted Hamiltonian where one component (beta) is rotated to a different basis in a selected ONV basis in which all configurations are selected (Full CI)
 *  Compare the evaluation of a direct matrix vector product to that of the matrix vector product evaluations and test the lowest eigenvalue against of the evaluated Hamiltonian a reference value
 */
BOOST_AUTO_TEST_CASE(Selected_H2O_Unrestricted) {

    // Psi4 and GAMESS' FCI energy (restricted)
    double reference_fci_energy = -75.0129803939602;

    // Create the molecular Hamiltonian in an AO basis
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::USpinorBasis<double, GQCP::GTOShell> spinor_basis {h2o, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spinor_basis, h2o);  // unrestricted Hamiltonian in the Löwdin basis

    // Transform the Hamiltonian to an orthonormal basis
    GQCP::basisTransform(spinor_basis, usq_hamiltonian, spinor_basis.lowdinOrthonormalizationMatrix().alpha());

    // Transform the beta component
    // Create stable unitairy matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes {usq_hamiltonian.spinHamiltonian(GQCP::Spin::alpha).core().parameters()};
    GQCP::basisTransform(spinor_basis, usq_hamiltonian, GQCP::TransformationMatrix<double>(saes.eigenvectors()), GQCP::Spin::beta);
    auto K = usq_hamiltonian.numberOfOrbitals() / 2;

    GQCP::SpinResolvedONVBasis fock_space {K, h2o.numberOfElectrons() / 2, h2o.numberOfElectrons() / 2};  // dim = 441
    GQCP::SpinResolvedSelectedONVBasis selected_fock_space {fock_space};

    // Evaluate the dense Hamiltonian
    GQCP::SquareMatrix<double> hamiltonian = selected_fock_space.evaluateOperatorDense(usq_hamiltonian, true);
    GQCP::SquareMatrix<double> hamiltonian_no_diagonal = selected_fock_space.evaluateOperatorDense(usq_hamiltonian, false);

    // Evaluate the diagonal
    GQCP::VectorX<double> hamiltonian_diagonal = selected_fock_space.evaluateOperatorDiagonal(usq_hamiltonian);

    // Evaluate the matvec
    GQCP::VectorX<double> matvec_evaluation = selected_fock_space.evaluateOperatorMatrixVectorProduct(usq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);

    // Evaluate the explicit matvec with the dense evaluations
    GQCP::VectorX<double> matvec_reference = hamiltonian * hamiltonian_diagonal;

    // Retrieve the lowest eigenvalue (FCI solution)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver {hamiltonian};

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
    double test_ci_energy = self_adjoint_eigensolver.eigenvalues()(0) + internuclear_repulsion_energy;

    // Test the energy with the reference
    BOOST_CHECK(std::abs(test_ci_energy - (reference_fci_energy)) < 1.0e-06);

    // Test if non-diagonal evaluation and diagonal evaluations are correct
    BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));

    // Test if the matvecs are identical
    BOOST_CHECK(matvec_evaluation.isApprox(matvec_reference));
}
