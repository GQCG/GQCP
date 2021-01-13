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

#define BOOST_TEST_MODULE "SpinResolvedONVBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/Transformations/transform.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


/**
 *  Test if the SpinResolvedONVBasis constructor throws when necessary.
 */
BOOST_AUTO_TEST_CASE(ProductONVBasis_constructor) {

    BOOST_CHECK_NO_THROW(GQCP::SpinResolvedONVBasis(10, 5, 5));
}


/**
 *  Check if the static SpinResolvedONVBasis basis dimension calculation is correct and if it can throw errors.
 */
BOOST_AUTO_TEST_CASE(ProductONVBasis_dimension) {

    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(10, 1, 1), 100);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(6, 2, 2), 225);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(8, 3, 3), 3136);

    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(10, 2, 0), 45);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(6, 3, 1), 120);
    BOOST_CHECK_EQUAL(GQCP::SpinResolvedONVBasis::calculateDimension(8, 4, 2), 1960);

    BOOST_CHECK_THROW(GQCP::SpinResolvedONVBasis::calculateDimension(60, 25, 25), std::overflow_error);
}


/**
 *  Check if the dense evaluation of a restricted one-electron operator matches between a full spin-resolved ONV basis and its 'selected' equivalent.
 * 
 *  The test system is H_6 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(evaluate_RSQOneElectronOperator_dense) {

    // Set up the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);


    // Set up the full spin-resolved ONV basis and its 'selected' equivalent.
    const auto K = sq_hamiltonian.numberOfOrbitals();
    const auto N_P = molecule.numberOfElectronPairs();
    GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};


    // Calculate the dense evaluations and check the result.
    const auto& h = sq_hamiltonian.core();
    const auto h_dense_specialized = onv_basis.evaluateOperatorDense(h);
    const auto h_dense_selected = selected_onv_basis.evaluateOperatorDense(h);

    BOOST_CHECK(h_dense_specialized.isApprox(h_dense_selected, 1.0e-12));
}


/**
 *  Check if the dense evaluation of a restricted two-electron operator matches between a full spin-resolved ONV basis and its 'selected' equivalent.
 * 
 *  The test system is a H_6 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(evaluate_RSQTwoElectronOperator_dense) {

    // Set up the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);


    // Set up the full spin-resolved ONV basis and its 'selected' equivalent.
    const auto K = sq_hamiltonian.numberOfOrbitals();
    const auto N_P = molecule.numberOfElectronPairs();
    GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};


    // Calculate the dense evaluations and check the result.
    const auto& g = sq_hamiltonian.twoElectron();
    const auto g_dense_specialized = onv_basis.evaluateOperatorDense(g);
    const auto g_dense_selected = selected_onv_basis.evaluateOperatorDense(g);

    BOOST_CHECK(g_dense_specialized.isApprox(g_dense_selected, 1.0e-12));
}


/**
 *  Check if the dense evaluation of a restricted Hamiltonian matches between a full spin-resolved ONV basis and its 'selected' equivalent.
 * 
 *  The test system is a H_6 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(evaluate_RSQHamiltonian_dense) {

    // Set up the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);

    // Set up the full spin-resolved ONV basis and its 'selected' equivalent.
    const auto K = hamiltonian.numberOfOrbitals();
    const auto N_P = molecule.numberOfElectronPairs();
    GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};


    // Calculate the dense evaluations and check the result.
    const auto H_dense_specialized = onv_basis.evaluateOperatorDense(hamiltonian);
    const auto H_dense_selected = selected_onv_basis.evaluateOperatorDense(hamiltonian);

    BOOST_CHECK(H_dense_specialized.isApprox(H_dense_selected, 1.0e-12));
}


/**
 *  Check if the diagonal of the matrix representation of a restricted one-electron operator is equal to the diagonal that is calculated through a specialized routine.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(restricted_one_electron_operator_diagonal) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-unresolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and the diagonal through a specialized routine, and check if they match.
    const auto dense_matrix = onv_basis.evaluateOperatorDense(hamiltonian.core());
    const auto diagonal_specialized = onv_basis.evaluateOperatorDiagonal(hamiltonian.core());

    BOOST_CHECK(diagonal_specialized.isApprox(dense_matrix.diagonal(), 1.0e-12));
}


/**
 *  Check if the diagonal of the matrix representation of a restricted Hamiltonian is equal to the diagonal that is calculated through a specialized routine.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(restricted_hamiltonian_diagonal) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-unresolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and the diagonal through a specialized routine, and check if they match.
    const auto dense_matrix = onv_basis.evaluateOperatorDense(hamiltonian);
    const auto diagonal_specialized = onv_basis.evaluateOperatorDiagonal(hamiltonian);

    BOOST_CHECK(diagonal_specialized.isApprox(dense_matrix.diagonal(), 1.0e-12));
}


/**
 *  Check if the matrix-vector product of a restricted one-electron operator through a direct evaluation (i.e. through the dense Hamiltonian matrix representation) and the specialized implementation are equal.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(restricted_one_electron_operator_dense_vs_matvec) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);

    // Set up the full spin-resolved ONV basis.
    const auto K = hamiltonian.numberOfOrbitals();
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and let it act on a random linear expansion.
    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis);
    const auto H_dense = onv_basis.evaluateOperatorDense(hamiltonian.core());
    const GQCP::VectorX<double> direct_mvp = H_dense * linear_expansion.coefficients();  // mvp: matrix-vector-product

    // Determine the specialized matrix-vector product and check if they are equal.
    const auto specialized_mvp = onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian.core(), linear_expansion.coefficients());

    BOOST_CHECK(specialized_mvp.isApprox(direct_mvp, 1.0e-08));
}


/**
 *  Check if the matrix-vector product of a restricted Hamiltonian through a direct evaluation (i.e. through the dense Hamiltonian matrix representation) and the specialized implementation are equal.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(restricted_Hamiltonian_dense_vs_matvec) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);

    // Set up the full spin-resolved ONV basis.
    const auto K = hamiltonian.numberOfOrbitals();
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and let it act on a random linear expansion.
    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis);
    const auto H_dense = onv_basis.evaluateOperatorDense(hamiltonian);
    const GQCP::VectorX<double> direct_mvp = H_dense * linear_expansion.coefficients();  // mvp: matrix-vector-product

    // Determine the specialized matrix-vector product and check if they are equal.
    const auto specialized_mvp = onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, linear_expansion.coefficients());

    BOOST_CHECK(specialized_mvp.isApprox(direct_mvp, 1.0e-08));
}


/**
 *  Check if the dense evaluation of an unrestricted Hamiltonian matches between a full spin-resolved ONV basis and its 'selected' equivalent.
 * 
 *  The test system is a H_6 in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(evaluate_USQHamiltonian_dense) {

    // Set up the molecular Hamiltonian in a random spin-orbital basis.
    const auto molecule = GQCP::Molecule::HChain(6, 0.742, 2);
    GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();
    hamiltonian.rotate(GQCP::UTransformation<double>::RandomUnitary(K));

    // Set up the full spin-resolved ONV basis and its 'selected' equivalent.
    const auto N_P = molecule.numberOfElectronPairs();
    GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};
    GQCP::SpinResolvedSelectedONVBasis selected_onv_basis {onv_basis};


    // Calculate the dense evaluations and check the result.
    const auto H_dense_specialized = onv_basis.evaluateOperatorDense(hamiltonian);
    const auto H_dense_selected = selected_onv_basis.evaluateOperatorDense(hamiltonian);

    BOOST_CHECK(H_dense_specialized.isApprox(H_dense_selected, 1.0e-12));
}


/**
 *  Check if the diagonal of the matrix representation of an unrestricted Hamiltonian is equal to the diagonal that is calculated through a specialized routine.
 * 
 *  The test system is a H6(2+)-chain with internuclear separation of 0.742 (a.u.) in an STO-3G basis.
 */
BOOST_AUTO_TEST_CASE(unrestricted_hamiltonian_diagonal) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::USpinorBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    const auto hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();

    // Set up the full spin-unresolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and the diagonal through a specialized routine, and check if they match.
    const auto dense_matrix = onv_basis.evaluateOperatorDense(hamiltonian);
    const auto diagonal_specialized = onv_basis.evaluateOperatorDiagonal(hamiltonian);

    BOOST_CHECK(diagonal_specialized.isApprox(dense_matrix.diagonal(), 1.0e-12));
}


/**
 *  Check if the matrix-vector product of an unrestricted Hamiltonian through a direct evaluation (i.e. through the dense Hamiltonian matrix representation) and the specialized implementation are equal.
 * 
 *  The test system is H2O in an STO-3G basisset, which has a FCI dimension of 441.
 */
BOOST_AUTO_TEST_CASE(unrestricted_dense_vs_matvec) {

    // Create the molecular Hamiltonian in a random unrestricted orthonormal spin-orbital basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();
    auto hamiltonian = GQCP::USQHamiltonian<double>::Molecular(spin_orbital_basis, molecule);
    const auto K = hamiltonian.numberOfOrbitals();
    hamiltonian.rotate(GQCP::UTransformation<double>::RandomUnitary(K));

    // Set up the full spin-resolved ONV basis.
    const GQCP::SpinResolvedONVBasis onv_basis {K, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Determine the Hamiltonian matrix and let it act on a random linear expansion.
    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::Random(onv_basis);
    const auto H_dense = onv_basis.evaluateOperatorDense(hamiltonian);
    const GQCP::VectorX<double> direct_mvp = H_dense * linear_expansion.coefficients();  // mvp: matrix-vector-product

    // Determine the specialized matrix-vector product and check if they are equal.
    const auto specialized_mvp = onv_basis.evaluateOperatorMatrixVectorProduct(hamiltonian, linear_expansion.coefficients());

    BOOST_CHECK(specialized_mvp.isApprox(direct_mvp, 1.0e-08));
}
