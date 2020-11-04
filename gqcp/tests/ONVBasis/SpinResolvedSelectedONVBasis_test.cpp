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
 *  Test the general functionality of the `expandWith` function, by testing throws and retrieving ONVs.
 */
BOOST_AUTO_TEST_CASE(expandWith) {

    // Set up an initial zero ONV basis.
    GQCP::SpinResolvedSelectedONVBasis onv_basis {3, 1, 1};


    // Add two compatible ONVs.
    BOOST_CHECK_NO_THROW(onv_basis.expandWith(SpinResolvedONV::FromString("001", "001")));
    BOOST_CHECK_NO_THROW(onv_basis.expandWith(SpinResolvedONV::FromString("010", "010")));


    // Check if we receive a throw if we add an ONV with an incompatible number of orbitals.
    BOOST_CHECK_THROW(onv_basis.expandWith(SpinResolvedONV::FromString("0001", "0100")), std::invalid_argument);


    // Check if we receive a throw if we add an ONV with an incompatible number of electrons.
    BOOST_CHECK_THROW(onv_basis.expandWith(SpinResolvedONV::FromString("011", "011")), std::invalid_argument);


    // Check if we can access the ONVs in order.
    BOOST_CHECK(onv_basis.onvWithIndex(0).asString() == "001|001");
    BOOST_CHECK(onv_basis.onvWithIndex(1).asString() == "010|010");
}


/**
 *  Evaluate the Hamiltonian in a selected ONV basis in which all configurations are selected (Full CI)
 *  Compare the evaluation of a direct matrix vector product to that of the matrix vector product evaluations and test the lowest eigenvalue against of the evaluated Hamiltonian a reference value
 */
// BOOST_AUTO_TEST_CASE(Selected_Evaluation_H2O) {

//     // Psi4 and GAMESS' FCI energy for H2O
//     double reference_fci_energy = -75.0129803939602;

//     // Create the molecular Hamiltonian in an AO basis
//     auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
//     GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {h2o, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in the Löwdin basis
//     auto K = sq_hamiltonian.numberOfOrbitals();

//     GQCP::SpinResolvedONVBasis onv_basis {K, h2o.numberOfElectrons() / 2, h2o.numberOfElectrons() / 2};  // dim = 441
//     GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};


//     // Evaluate the dense Hamiltonian
//     GQCP::SquareMatrix<double> hamiltonian = selected_onv_basis.evaluateOperatorDense(sq_hamiltonian, true);
//     GQCP::SquareMatrix<double> hamiltonian_no_diagonal = selected_onv_basis.evaluateOperatorDense(sq_hamiltonian, false);

//     // Evaluate the diagonal
//     GQCP::VectorX<double> hamiltonian_diagonal = selected_onv_basis.evaluateOperatorDiagonal(sq_hamiltonian);

//     // Evaluate the matvec
//     GQCP::VectorX<double> matvec_evaluation = selected_onv_basis.evaluateOperatorMatrixVectorProduct(sq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);

//     // Calculate the explicit matvec with the dense evaluations
//     GQCP::VectorX<double> matvec_reference = hamiltonian * hamiltonian_diagonal;

//     // Retrieve the lowest eigenvalue (FCI solution)
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver {hamiltonian};

//     double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
//     double test_energy = self_adjoint_eigensolver.eigenvalues()(0) + internuclear_repulsion_energy;

//     // Test the energy with the reference
//     BOOST_CHECK(std::abs(test_energy - reference_fci_energy) < 1e-6);

//     // Test if the non-diagonal evaluation and diagonal evaluations are correct
//     BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));

//     // Test if the matvecs are identical
//     BOOST_CHECK(matvec_evaluation.isApprox(matvec_reference));
// }


/**
 *  Evaluate the an unrestricted Hamiltonian where one component (beta) is rotated to a different basis in a selected ONV basis in which all configurations are selected (Full CI)
 *  Compare the evaluation of a direct matrix vector product to that of the matrix vector product evaluations and test the lowest eigenvalue against of the evaluated Hamiltonian a reference value
 */
// BOOST_AUTO_TEST_CASE(Selected_H2O_Unrestricted) {

//     // Psi4 and GAMESS' FCI energy (restricted)
//     double reference_fci_energy = -75.0129803939602;

//     // Create the molecular Hamiltonian in an AO basis
//     auto h2o = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
//     GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {h2o, "STO-3G"};
//     spinor_basis.lowdinOrthonormalize();
//     auto usq_hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spinor_basis, h2o);  // unrestricted Hamiltonian in the Löwdin basis

//     // Transform the Hamiltonian to an orthonormal basis
//     GQCP::basisTransform(spinor_basis, usq_hamiltonian, spinor_basis.lowdinOrthonormalizationMatrix().alpha());

//     // Transform the beta component
//     // Create stable unitairy matrix
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes {usq_hamiltonian.spinHamiltonian(GQCP::Spin::alpha).core().parameters()};
//     GQCP::basisTransform(spinor_basis, usq_hamiltonian, GQCP::TransformationMatrix<double>(saes.eigenvectors()), GQCP::Spin::beta);
//     auto K = usq_hamiltonian.numberOfOrbitals() / 2;

//     GQCP::SpinResolvedONVBasis onv_basis {K, h2o.numberOfElectrons() / 2, h2o.numberOfElectrons() / 2};  // dim = 441
//     GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};

//     // Evaluate the dense Hamiltonian
//     GQCP::SquareMatrix<double> hamiltonian = selected_onv_basis.evaluateOperatorDense(usq_hamiltonian, true);
//     GQCP::SquareMatrix<double> hamiltonian_no_diagonal = selected_onv_basis.evaluateOperatorDense(usq_hamiltonian, false);

//     // Evaluate the diagonal
//     GQCP::VectorX<double> hamiltonian_diagonal = selected_onv_basis.evaluateOperatorDiagonal(usq_hamiltonian);

//     // Evaluate the matvec
//     GQCP::VectorX<double> matvec_evaluation = selected_onv_basis.evaluateOperatorMatrixVectorProduct(usq_hamiltonian, hamiltonian_diagonal, hamiltonian_diagonal);

//     // Evaluate the explicit matvec with the dense evaluations
//     GQCP::VectorX<double> matvec_reference = hamiltonian * hamiltonian_diagonal;

//     // Retrieve the lowest eigenvalue (FCI solution)
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> self_adjoint_eigensolver {hamiltonian};

//     // Calculate the total FCI energy
//     double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2o).value();
//     double test_ci_energy = self_adjoint_eigensolver.eigenvalues()(0) + internuclear_repulsion_energy;

//     // Test the energy with the reference
//     BOOST_CHECK(std::abs(test_ci_energy - (reference_fci_energy)) < 1.0e-06);

//     // Test if non-diagonal evaluation and diagonal evaluations are correct
//     BOOST_CHECK(hamiltonian.isApprox(hamiltonian_no_diagonal + GQCP::SquareMatrix<double>(hamiltonian_diagonal.asDiagonal())));

//     // Test if the matvecs are identical
//     BOOST_CHECK(matvec_evaluation.isApprox(matvec_reference));
// }
