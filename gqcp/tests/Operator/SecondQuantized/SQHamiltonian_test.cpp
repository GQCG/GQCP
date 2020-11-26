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

#define BOOST_TEST_MODULE "HamiltonianParameters"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/miscellaneous.hpp"

#include <boost/math/constants/constants.hpp>


/*
 *  HELPER FUNCTIONS
 */

/**
 *  @return a toy 2-DM:
 *      d(p,q,r,s) = delta_pq delta_rs
 */
GQCP::Orbital2DM<double> calculateToy2DM() {

    auto d = GQCP::Orbital2DM<double>::Zero(2);

    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if ((p == q) && (r == s)) {
                        d(p, q, r, s) = 1;
                    }
                }
            }
        }
    }

    return d;
};


/**
 *  @return toy 2-electron integrals:
 *      g(p,q,r,s) = -0.5 delta_pq delta_rs
 */
GQCP::ScalarRSQTwoElectronOperator<double> calculateToyTwoElectronIntegrals() {

    GQCP::SquareRankFourTensor<double> g {2};
    g.setZero();

    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if ((p == q) && (r == s)) {
                        g(p, q, r, s) = -0.5;
                    }
                }
            }
        }
    }

    return GQCP::ScalarRSQTwoElectronOperator<double> {g};
};


/*
 *  UNIT TESTS - CONSTRUCTORS
 */


/**
 *  Check if the constructor works as expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {

    const size_t K = 3;

    // Create one- and two-electron operators with compatible dimensions.
    const auto h_core = GQCP::ScalarRSQOneElectronOperator<double>::Random(K);
    const auto g = GQCP::ScalarRSQTwoElectronOperator<double>::Random(K);


    // Check if a correct constructor works.
    BOOST_CHECK_NO_THROW(GQCP::RSQHamiltonian<double> sq_hamiltonian(h_core, g));


    // Check if wrong arguments result in a throw.
    const auto h_core_faulty = GQCP::ScalarRSQOneElectronOperator<double>::Random(K + 1);
    const auto g_faulty = GQCP::ScalarRSQTwoElectronOperator<double>::Random(K + 1);

    BOOST_CHECK_THROW(GQCP::RSQHamiltonian<double> sq_hamiltonian(h_core, g_faulty), std::invalid_argument);
    BOOST_CHECK_THROW(GQCP::RSQHamiltonian<double> sq_hamiltonian(h_core_faulty, g), std::invalid_argument);

    BOOST_CHECK_NO_THROW(GQCP::RSQHamiltonian<double> sq_hamiltonian(h_core_faulty, g_faulty));
}


/*
 *  UNIT TESTS - NAMED CONSTRUCTORS
 */

BOOST_AUTO_TEST_CASE(constructMolecularHamiltonianParameters) {

    // Set up a basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_szabo.xyz");


    // Check if we can construct the molecular Hamiltonian
    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {h2, "STO-3G"};
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto g = sq_hamiltonian.twoElectron().parameters();


    // Check with reference values from Szabo
    GQCP::SquareMatrix<double> ref_S {2};
    // clang-format off
    ref_S << 1.0,    0.6593,
             0.6593, 1.0;
    // clang-format on

    GQCP::SquareMatrix<double> ref_H_core {2};
    // clang-format off
    ref_H_core << -1.1204, -0.9584,
                  -0.9584, -1.1204;
    // clang-format on

    BOOST_CHECK(spinor_basis.overlap().parameters().isApprox(ref_S, 1.0e-04));
    BOOST_CHECK(sq_hamiltonian.core().parameters().isApprox(ref_H_core, 1.0e-04));

    BOOST_CHECK(std::abs(g(0, 0, 0, 0) - 0.7746) < 1.0e-04);
    BOOST_CHECK(std::abs(g(0, 0, 0, 0) - g(1, 1, 1, 1)) < 1.0e-12);
    BOOST_CHECK(std::abs(g(0, 0, 1, 1) - 0.5697) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1, 0, 0, 0) - 0.4441) < 1.0e-04);
    BOOST_CHECK(std::abs(g(1, 0, 0, 0) - g(1, 1, 1, 0)) < 1.0e-12);
    BOOST_CHECK(std::abs(g(1, 0, 1, 0) - 0.2970) < 1.0e-04);
}


BOOST_AUTO_TEST_CASE(FCIDUMP_reader) {

    auto fcidump_ham_par = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/beh_cation_631g_caitlin.FCIDUMP");

    // Check if the one-electron integrals are read in correctly from a previous implementation
    GQCP::SquareMatrix<double> h_SO = fcidump_ham_par.core().parameters();

    BOOST_CHECK(std::abs(h_SO(0, 0) - (-8.34082)) < 1.0e-5);
    BOOST_CHECK(std::abs(h_SO(5, 1) - 0.381418) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(14, 0) - 0.163205) < 1.0e-6);
    BOOST_CHECK(std::abs(h_SO(13, 6) - (-5.53204e-16)) < 1.0e-16);
    BOOST_CHECK(std::abs(h_SO(15, 11) - (-0.110721)) < 1.0e-6);


    // Check if the two-electron integrals are read in correctly from a previous implementation
    GQCP::SquareRankFourTensor<double> g_SO = fcidump_ham_par.twoElectron().parameters();

    BOOST_CHECK(std::abs(g_SO(2, 5, 4, 4) - 0.0139645) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(2, 6, 3, 0) - 5.16622e-18) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(3, 1, 3, 0) - (-0.0141251)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4, 6, 4, 6) - 0.0107791) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(4, 15, 11, 1) - (9.33375e-19)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(6, 10, 5, 9) - (-3.81422e-18)) < 1.0e-17);
    BOOST_CHECK(std::abs(g_SO(7, 7, 2, 1) - (-0.031278)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(8, 15, 9, 9) - (-2.80093e-17)) < 1.0e-16);
    BOOST_CHECK(std::abs(g_SO(9, 14, 0, 9) - 0.00161985) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(10, 1, 4, 3) - 0.00264603) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(11, 4, 9, 3) - (-0.0256623)) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(12, 9, 0, 4) - 0.0055472) < 1.0e-6);
    BOOST_CHECK(std::abs(g_SO(13, 15, 15, 13) - 0.00766898) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(14, 2, 12, 3) - 0.0104266) < 1.0e-7);
    BOOST_CHECK(std::abs(g_SO(15, 5, 10, 10) - 0.00562608) < 1.0e-7);
}


BOOST_AUTO_TEST_CASE(FCIDUMP_reader_HORTON) {

    // Check the same reference value that HORTON does
    auto fcidump_ham_par = GQCP::RSQHamiltonian<double>::FromFCIDUMP("data/h2_psi4_horton.FCIDUMP");

    GQCP::SquareRankFourTensor<double> g_SO = fcidump_ham_par.twoElectron().parameters();
    BOOST_CHECK(std::abs(g_SO(6, 5, 1, 0) - 0.0533584656) < 1.0e-7);
}


/*
 *  UNIT TESTS - METHODS
 */

BOOST_AUTO_TEST_CASE(calculate_generalized_Fock_matrix_and_super_invalid_arguments) {

    // Initialize toy HamiltonianParameters
    GQCP::SquareMatrix<double> h = GQCP::SquareMatrix<double>::Zero(2);
    GQCP::SquareRankFourTensor<double> g {2};
    GQCP::RSQHamiltonian<double> sq_hamiltonian {GQCP::ScalarRSQOneElectronOperator<double>(h), GQCP::ScalarRSQTwoElectronOperator<double>(g)};


    // Create valid and invalid density matrices (with respect to the dimensions of the SOBasis)
    GQCP::Orbital1DM<double> D_valid = GQCP::Orbital1DM<double>::Zero(2);
    GQCP::Orbital1DM<double> D_invalid = GQCP::Orbital1DM<double>::Zero(3);

    GQCP::Orbital2DM<double> d_valid {2};
    GQCP::Orbital2DM<double> d_invalid {3};


    // Test a faulty function calls
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateFockianMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateFockianMatrix(D_valid, d_invalid), std::invalid_argument);

    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateSuperFockianMatrix(D_invalid, d_valid), std::invalid_argument);
    BOOST_REQUIRE_THROW(sq_hamiltonian.calculateSuperFockianMatrix(D_valid, d_invalid), std::invalid_argument);


    // Test correct function calls
    sq_hamiltonian.calculateFockianMatrix(D_valid, d_valid);
    sq_hamiltonian.calculateSuperFockianMatrix(D_valid, d_valid);
}


BOOST_AUTO_TEST_CASE(calculate_Fockian_and_super) {

    // We test the function by a manual calculation of nonsensical toy 1- and 2-DMs and one- and two-electron integrals
    // Set up the toy 1- and 2-DMs
    GQCP::Orbital1DM<double> D {2};
    // clang-format off
    D << 0, 1,
         1, 0;
    // clang-format on

    auto d = calculateToy2DM();

    // Set up the toy Hamiltonian
    GQCP::SquareMatrix<double> h = GQCP::SquareMatrix<double>::Zero(2);
    // clang-format off
    h << 1, 0,
         0, 1;
    // clang-format on

    auto g = calculateToyTwoElectronIntegrals();
    GQCP::RSQHamiltonian<double> sq_hamiltonian {GQCP::ScalarRSQOneElectronOperator<double>(h), g};


    // Construct the reference Fockian matrix
    GQCP::SquareMatrix<double> F_ref {2};
    // clang-format off
    F_ref << -1.00,  1.00,
              1.00, -1.00;
    // clang-format on

    // Construct the reference super generalized Fock matrix
    GQCP::SquareRankFourTensor<double> G_ref {2};
    G_ref.setZero();
    for (size_t p = 0; p < 2; p++) {
        for (size_t q = 0; q < 2; q++) {
            for (size_t r = 0; r < 2; r++) {
                for (size_t s = 0; s < 2; s++) {
                    if (q == r) {
                        G_ref(p, q, r, s) += F_ref(p, s);
                    }

                    // One-electron part is simplified by manual calculation
                    if (p == s) {
                        G_ref(p, q, r, s) -= D(r, q);
                    }

                    // Two-electron part is simplified by manual calculation
                    if ((p == s) && (q == r)) {
                        G_ref(p, q, r, s) += 1;
                    }
                }
            }
        }
    }


    BOOST_CHECK(F_ref.isApprox(sq_hamiltonian.calculateFockianMatrix(D, d), 1.0e-12));
    BOOST_CHECK(G_ref.isApprox(sq_hamiltonian.calculateSuperFockianMatrix(D, d), 1.0e-12));
}


/**
 *  Check if the implementation of the Edmiston-Ruedenberg localization index works as expected.
 */
BOOST_AUTO_TEST_CASE(calculateEdmistonRuedenbergLocalizationIndex) {

    // Create a toy Hamiltonian: only the two-electron integrals are important
    const size_t K = 5;
    const GQCP::SquareMatrix<double> H_op = GQCP::SquareMatrix<double>::Random(K);

    GQCP::SquareRankFourTensor<double> g_op {K};
    g_op.setZero();
    for (size_t p = 0; p < K; p++) {
        g_op(p, p, p, p) = 2 * static_cast<float>(p);
    }

    GQCP::RSQHamiltonian<double> sq_hamiltonian {GQCP::ScalarRSQOneElectronOperator<double>(H_op), GQCP::ScalarRSQTwoElectronOperator<double>(g_op)};


    // Check the values for the Edmiston-Ruedenberg localization index.
    const auto orbital_space1 = GQCP::OrbitalSpace::Implicit({{GQCP::OccupationType::k_occupied, 3}});  // 3 occupied spatial orbitals
    const auto orbital_space2 = GQCP::OrbitalSpace::Implicit({{GQCP::OccupationType::k_occupied, 4}});  // 3 occupied spatial orbitals

    BOOST_CHECK(std::abs(sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(orbital_space1) - 6.0) < 1.0e-08);
    BOOST_CHECK(std::abs(sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(orbital_space2) - 12.0) < 1.0e-08);
}


/**
 *  Test if we can succesfully initialize NO+ at long intra molecular distance.
 */
BOOST_AUTO_TEST_CASE(dissociatedMoleculeParameters) {

    GQCP::Nucleus N {7, 3.5, 0, 0};
    GQCP::Nucleus O {8, -3.5, 0, 0};
    std::vector<GQCP::Nucleus> nuclei {N, O};
    GQCP::Molecule NO {nuclei, +1};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {NO, "STO-3G"};
    BOOST_CHECK_NO_THROW(GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, NO));
}
