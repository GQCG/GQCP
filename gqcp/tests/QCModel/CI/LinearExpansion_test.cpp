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

#define BOOST_TEST_MODULE "LinearExpansion"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCModel/CI/LinearExpansion.hpp"


/**
 *  Test if a GAMESS-US expansion file is correctly read in.
 */
BOOST_AUTO_TEST_CASE(reader_test) {

    // Provide the reference values.
    const GQCP::VectorX<double> ref_coefficients = GQCP::VectorX<double>::Unit(2, 0);  // (size, position)
    const std::string alpha1_ref = "0000000000000000000000000000000000000000000001";
    const std::string alpha2_ref = "0000000000000000000000000000000000000000000001";
    const std::string beta1_ref = "0000000000000000000000000000000000000000000001";
    const std::string beta2_ref = "0000000000000000000000000000000000000000000010";


    // Read in the GAMESS-US file and check the results.
    const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedSelectedONVBasis>::FromGAMESSUS("data/test_GAMESS_expansion");

    // Check if the expansion coefficients are correct.
    BOOST_CHECK(linear_expansion.coefficients().isApprox(ref_coefficients, 1.0e-08));

    // Check if the parsed ONVs are correct.
    const auto onv1 = linear_expansion.onvBasis().onvWithIndex(0);
    const auto onv1_alpha = onv1.onv(GQCP::Spin::alpha).asString();
    const auto onv1_beta = onv1.onv(GQCP::Spin::beta).asString();

    BOOST_CHECK(onv1_alpha == alpha1_ref);
    BOOST_CHECK(onv1_beta == beta1_ref);

    const auto onv2 = linear_expansion.onvBasis().onvWithIndex(1);
    const auto onv2_alpha = onv2.onv(GQCP::Spin::alpha).asString();
    const auto onv2_beta = onv2.onv(GQCP::Spin::beta).asString();

    BOOST_CHECK(onv2_alpha == alpha2_ref);
    BOOST_CHECK(onv2_beta == beta2_ref);
}


/**
 *  Check if the calculation of the Shannon entropy is correctly implemented by comparing with a manual calculation.
 */
BOOST_AUTO_TEST_CASE(shannon_entropy) {

    // Set up a test spin-resolved ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis {8, 3};  // 8 spinors, 3 electrons


    // Check the Shannon entropy of a Hartree-Fock expansion
    const auto hartree_fock_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::HartreeFock(onv_basis);
    BOOST_CHECK(hartree_fock_expansion.calculateShannonEntropy() < 1.0e-12);  // should be 0


    // Check the maximal entropy, corresponding to a wave function with all equal coefficients different from zero.
    const auto constant_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Constant(onv_basis);
    const double reference_entropy = std::log2(onv_basis.dimension());  // manual derivation
    BOOST_CHECK(std::abs(constant_expansion.calculateShannonEntropy() - reference_entropy) < 1.0e-12);
}


/**
 *  Check if the basis transformation of a linear expansion inside the full spin-resolved ONV basis is correctly implemented: we compare the direct transformation of the expansion coefficients with another FCI calculation using the transformed spinor basis.
 *  The test system is a linear H chain H3-//STO-3G, with an internuclear charge 0.742 bohr.
 */
BOOST_AUTO_TEST_CASE(transform_wave_function_h3) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(3, 0.742, -1);  // charge -1
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation and calculate the transformation of the linear expansion coefficients.
    const auto U_random = GQCP::RTransformation<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransform(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::rotate(U_random, spinor_basis, sq_hamiltonian);

    auto environment_indirect = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_indirect = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_indirect = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_indirect, environment_indirect).groundStateParameters();
    BOOST_CHECK(linear_expansion_direct.isApprox(linear_expansion_indirect, 1.0e-12));
}


/**
 *  Check if the basis transformation of a linear expansion inside the full spin-resolved ONV basis is correctly implemented: we compare the direct transformation of the expansion coefficients with another FCI calculation using the transformed spinor basis.
 *  The test system is a linear H chain H4//STO-3G, with an internuclear charge 0.742 bohr.
 */
BOOST_AUTO_TEST_CASE(transform_wave_function_h4) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(4, 0.742);
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation and calculate the transformation of the linear expansion coefficients.
    const auto U_random = GQCP::RTransformation<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransform(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::rotate(U_random, spinor_basis, sq_hamiltonian);

    auto environment_indirect = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_indirect = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_indirect = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_indirect, environment_indirect).groundStateParameters();
    BOOST_CHECK(linear_expansion_direct.isApprox(linear_expansion_indirect, 1.0e-12));
}


/**
 *  Check if the basis transformation of a linear expansion inside the full spin-resolved ONV basis is correctly implemented: we compare the direct transformation of the expansion coefficients with another FCI calculation using the transformed spinor basis.
 *  The test system is a linear H chain H5//STO-3G, with an internuclear charge 0.742 bohr.
 */
BOOST_AUTO_TEST_CASE(transform_wave_function_h5) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(5, 0.742);
    const auto N_alpha = 3;
    const auto N_beta = 2;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation and calculate the transformation of the linear expansion coefficients.
    const auto U_random = GQCP::RTransformation<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransform(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::rotate(U_random, spinor_basis, sq_hamiltonian);

    auto environment_indirect = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_indirect = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_indirect = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_indirect, environment_indirect).groundStateParameters();
    BOOST_CHECK(linear_expansion_direct.isApprox(linear_expansion_indirect, 1.0e-12));
}


/**
 *  Test if the LinearExpansions generated by a SpinUnresolvedONVBasis basis are normalized.
 */
BOOST_AUTO_TEST_CASE(expansions) {

    // Create a toy spin-unresolved ONV basis.
    GQCP::SpinUnresolvedONVBasis onv_basis {8, 3};


    // Check if ::Constant yields a normalized linear expansion.
    const auto constant_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Constant(onv_basis);
    BOOST_CHECK(std::abs(constant_expansion.coefficients().norm() - 1.0) < 1.0e-12);  // check if the coefficient vector is normalized


    // Check if ::HartreeFock yields a normalized linear expansion.
    const auto hartree_fock_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::HartreeFock(onv_basis);
    BOOST_CHECK(std::abs(hartree_fock_expansion.coefficients().norm() - 1.0) < 1.0e-12);  // check if the coefficient vector is normalized
    BOOST_CHECK(std::abs(hartree_fock_expansion.coefficients()(0) - 1.0) < 1.0e-12);      // the Hartree-Fock determinant should be the first one


    // Check if ::Random yields a normalized linear expansion.
    const auto random_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis);
    BOOST_CHECK(std::abs(random_expansion.coefficients().norm() - 1.0) < 1.0e-12);  // check if the coefficient vector is normalized
}


/**
 *  Check if calculateNDMElement throws as expected.
 */
BOOST_AUTO_TEST_CASE(calculateNDMElement_throws) {

    // Set up an example linear expansion.
    const size_t M = 3;
    const size_t N = 1;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    GQCP::VectorX<double> coefficients {onv_basis.dimension()};
    coefficients << 1, 2, -3;

    const GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> linear_expansion {onv_basis, coefficients};


    // Check if calculateNDMElement throws as expected.
    BOOST_CHECK_THROW(linear_expansion.calculateNDMElement({3}, {0}), std::invalid_argument);  // bra-index is out of bounds
    BOOST_CHECK_THROW(linear_expansion.calculateNDMElement({0}, {3}), std::invalid_argument);  // ket-index is out of bounds
}


/**
 *  Check some 1-DM values calculated through the general function calculateNDMElement.
 */
BOOST_AUTO_TEST_CASE(calculateNDMElement_1DM) {

    // Set up an example linear expansion.
    const size_t M = 3;
    const size_t N = 1;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    GQCP::VectorX<double> coefficients {onv_basis.dimension()};
    coefficients << 1, 2, -3;

    const GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> linear_expansion {onv_basis, coefficients};


    // Check some 1-DM values.
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0}, {0}) - 1.0) < 1.0e-12);     // d(0,0) : a^\dagger_0 a_0
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0}, {1}) - 2.0) < 1.0e-12);     // d(0,1) : a^\dagger_0 a_1
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({2}, {1}) - (-6.0)) < 1.0e-12);  // d(2,1) : a^\dagger_2 a_1
}


/**
 *  Check some 2-DM values calculated through the general function calculateNDMElement.
 */
BOOST_AUTO_TEST_CASE(calculateNDMElement_2DM) {

    // Set up an example linear expansion.
    const size_t M = 3;
    const size_t N = 2;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    GQCP::VectorX<double> coefficients {onv_basis.dimension()};
    coefficients << 1, 2, -3;

    const GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> linear_expansion {onv_basis, coefficients};


    // Check some 2-DM values.
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 1}, {2, 1}) - (-3.0)) < 1.0e-12);  // d(0,1,1,2) : a^\dagger_0 a^\dagger_1 a_2 a_1
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({2, 0}, {1, 0}) - (-2.0)) < 1.0e-12);  // d(2,0,0,1) : a^\dagger_2 a^\dagger_0 a^1 a_0
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 2}, {0, 2}) - (-4.0)) < 1.0e-12);  // d(0,2,2,0) : a^\dagger_0 a^dagger_2 a_0 a_2
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 0}, {0, 2}) - 0.0) < 1.0e-12);     // d(0,2,0,0) : a^\dagger_0 a^dagger_0 a_0 a_2, double annihilation gives 0.0
}


/**
 *  Check some 3-DM values calculated through the general function calculateNDMElement.
 */
BOOST_AUTO_TEST_CASE(calculateNDMElement_3DM) {

    // Set up an example linear expansion.
    const size_t M = 5;
    const size_t N = 4;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    GQCP::VectorX<double> coefficients {onv_basis.dimension()};
    coefficients << 1, 1, -2, 4, -5;

    const GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> linear_expansion {onv_basis, coefficients};


    // Check some 3-DM values.
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 0, 1}, {1, 0, 2}) - 0.0) < 1.0e-12);  // zero because two times the same index
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({1, 0, 3}, {4, 1, 2}) - 0.0) < 1.0e-12);  // zero because no fully annihilated bras and kets match
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 1, 2}, {2, 1, 0}) - 2.0) < 1.0e-12);
    BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0, 1, 2}, {0, 1, 3}) - 2.0) < 1.0e-12);
}


// /**
//  *  Check if the projection of |UHF> and |GHF> (the equivalent ONV in a generalized spinor basis) onto |RHF> is the same.
//  */
// BOOST_AUTO_TEST_CASE(overlap_GHF_UHF_on_RHF) {

//     const auto molecule = GQCP::Molecule::HRingFromDistance(4, 1.0);

//     GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
//     GQCP::SquareMatrix<double> C {4};  // RHF canonical orbitals for this system (Xeno)
//     // clang-format off
//     C << -0.27745359, -0.8505133,   0.85051937,  2.02075317,
//          -0.27745362, -0.85051937, -0.8505133,  -2.02075317,
//          -0.27745359,  0.8505133,  -0.85051937,  2.02075317,
//          -0.27745362,  0.85051937,  0.8505133,  -2.02075317;
//     // clang-format on
//     r_spinor_basis.transform(C);


//     // GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> u_spinor_basis {molecule, "STO-3G"};
//     GQCP::SquareMatrix<double> C_alpha {4};  // UHF alpha canonical orbitals for this system (Xeno), triplet
//     // clang-format off
//     C_alpha << -1.75646828e-01, -1.20606646e-06,  1.20281173e+00,  2.03213486e+00,
//                -3.78560533e-01, -1.20281173e+00, -1.20606647e-06, -2.00427438e+00,
//                -1.75646828e-01,  1.20606646e-06, -1.20281173e+00,  2.03213486e+00,
//                -3.78560533e-01,  1.20281173e+00,  1.20606646e-06, -2.00427438e+00;
//     // clang-format on

//     GQCP::SquareMatrix<double> C_beta {4};  // UHF alpha canonical orbitals for this system (Xeno), triplet
//     // clang-format off
//     C_beta << -3.78560533e-01,  1.20281173e+00,  1.21724557e-06,  2.00427438e+00,
//               -1.75646828e-01,  1.21724558e-06, -1.20281173e+00, -2.03213486e+00,
//               -3.78560533e-01, -1.20281173e+00, -1.21724558e-06,  2.00427438e+00,
//               -1.75646828e-01, -1.21724558e-06,  1.20281173e+00, -2.03213486e+00;
//     // clang-format on
//     // u_spinor_basis.transform(C_alpha, C_beta);


//     const auto u_spinor_basis = GQCP::USpinOrbitalBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);


//     const auto uhf_onv = GQCP::SpinResolvedONV::UHF(4, 2, 2);
//     const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::FromONVProjection(uhf_onv, r_spinor_basis, u_spinor_basis);

//     std::cout << linear_expansion.coefficients() << std::endl;
// }


/**
 *  Check if the 1- and 2-DMs for a full spin-resolved ONV basis are equal to the 'selected' case.
 *
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 441. However, we're choosing a different number of alpha and beta electrons. (N_alpha = 4, N_beta = 6)
 */
BOOST_AUTO_TEST_CASE(spin_resolved_vs_spin_resolved_selected_DMs) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const size_t N_alpha = 4;
    const size_t N_beta = 6;

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_specialized = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Convert the 'specialized' linear expansion into a 'selected' linear expansion.
    const GQCP::SpinResolvedSelectedONVBasis onv_basis_selected {onv_basis};
    const auto linear_expansion_selected = GQCP::LinearExpansion<GQCP::SpinResolvedSelectedONVBasis>(onv_basis_selected, linear_expansion_specialized.coefficients());


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto D_specialized = linear_expansion_specialized.calculateSpinResolved1DM();
    const auto D_selected = linear_expansion_selected.calculateSpinResolved1DM();

    BOOST_CHECK(D_specialized.orbitalDensity().isApprox(D_selected.orbitalDensity(), 1.0e-12));
    BOOST_CHECK(D_specialized.alpha().isApprox(D_selected.alpha(), 1.0e-12));
    BOOST_CHECK(D_specialized.beta().isApprox(D_selected.beta(), 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto d_specialized = linear_expansion_specialized.calculateSpinResolved2DM();
    const auto d_selected = linear_expansion_selected.calculateSpinResolved2DM();

    BOOST_CHECK(d_specialized.alphaAlpha().isApprox(d_selected.alphaAlpha(), 1.0e-12));
    BOOST_CHECK(d_specialized.alphaBeta().isApprox(d_selected.alphaBeta(), 1.0e-12));
    BOOST_CHECK(d_specialized.betaAlpha().isApprox(d_selected.betaAlpha(), 1.0e-12));
    BOOST_CHECK(d_specialized.betaBeta().isApprox(d_selected.betaBeta(), 1.0e-12));
    BOOST_CHECK(d_specialized.orbitalDensity().isApprox(d_selected.orbitalDensity(), 1.0e-12));
}


/**
 *  Check if the 1- and 2-DMs for a seniority-zero ONV basis are equal to the 'selected' case.
 *
 *  The system of interested is H2O//STO-3G, with 7 spatial orbitals and a Fock space dimension of 21.
 */
BOOST_AUTO_TEST_CASE(seniority_zero_vs_spin_resolved_selected_DMs) {

    // Set up the molecular Hamiltonian in a Löwdin-orthonormalized spinor basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_Psi4_GAMESS.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();
    spinor_basis.lowdinOrthonormalize();

    const auto sq_hamiltonian = GQCP::RSQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense DOCI calculation.
    const GQCP::SeniorityZeroONVBasis onv_basis {K, N_P};

    auto environment = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_specialized = GQCP::QCMethod::CI<GQCP::SeniorityZeroONVBasis>(onv_basis).optimize(solver, environment).groundStateParameters();


    // Convert the 'specialized' linear expansion into a 'selected' linear expansion.
    const GQCP::SpinResolvedSelectedONVBasis onv_basis_selected {onv_basis};
    const auto linear_expansion_selected = GQCP::LinearExpansion<GQCP::SpinResolvedSelectedONVBasis>(onv_basis_selected, linear_expansion_specialized.coefficients());


    // Calculate the 1-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto D_specialized = linear_expansion_specialized.calculateSpinResolved1DM();
    const auto D_selected = linear_expansion_selected.calculateSpinResolved1DM();

    BOOST_CHECK(D_specialized.orbitalDensity().isApprox(D_selected.orbitalDensity(), 1.0e-12));
    BOOST_CHECK(D_specialized.alpha().isApprox(D_selected.alpha(), 1.0e-12));
    BOOST_CHECK(D_specialized.beta().isApprox(D_selected.beta(), 1.0e-12));


    // Calculate the 2-DMs using specialized spin-resolved and 'selected' routines, and check if they are equal.
    const auto d_specialized = linear_expansion_specialized.calculateSpinResolved2DM();
    const auto d_selected = linear_expansion_selected.calculateSpinResolved2DM();

    BOOST_CHECK(d_specialized.alphaAlpha().isApprox(d_selected.alphaAlpha(), 1.0e-12));
    BOOST_CHECK(d_specialized.alphaBeta().isApprox(d_selected.alphaBeta(), 1.0e-12));
    BOOST_CHECK(d_specialized.betaAlpha().isApprox(d_selected.betaAlpha(), 1.0e-12));
    BOOST_CHECK(d_specialized.betaBeta().isApprox(d_selected.betaBeta(), 1.0e-12));
    BOOST_CHECK(d_specialized.orbitalDensity().isApprox(d_selected.orbitalDensity(), 1.0e-12));
}


// /**
//  *  Check some 1-DM values calculated for the SpinUnresolvedONVBasis by comparing them to the general function calculateNDMElement.
//  */
// BOOST_AUTO_TEST_CASE(calculate1DM_SpinUnresolved_NDM) {

//     // Set up an example linear expansion in a spin-unresolved ONV basis.
//     const size_t M = 5;
//     const size_t N = 2;
//     const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

//     const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis);

//     // Check some 1-DM values with calculateNDMElement().
//     const auto D_specialized = linear_expansion.calculate1DM();
//     BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0}, {0}) - D_specialized(0, 0)) < 1.0e-12);  // D(0,0) : a^\dagger_0 a_0
//     BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({0}, {1}) - D_specialized(0, 1)) < 1.0e-12);  // D(0,1) : a^\dagger_0 a_1
//     BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({2}, {1}) - D_specialized(2, 1)) < 1.0e-12);  // D(2,1) : a^\dagger_2 a_1
//     BOOST_CHECK(std::abs(linear_expansion.calculateNDMElement({4}, {4}) - D_specialized(4, 4)) < 1.0e-12);  // D(4,4) : a^\dagger_4 a_4
// }


/**
 *  Check some 1-DM values calculated for the SpinUnresolvedONVBasis by comparing them to an equivalent spin-resolved calculation.
 */
BOOST_AUTO_TEST_CASE(calculate1DM_SpinUnresolved_ref_SpinResolved) {

    // Set up an example linear expansion in a spin-unresolved ONV basis.
    const size_t M = 5;
    const size_t N = 2;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    const auto linear_expansion_unresolved = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis);
    const auto D_unresolved = linear_expansion_unresolved.calculate1DM();


    // Create an equivalent, spin-resolved linear expansion.
    const GQCP::SpinResolvedONVBasis onv_basis_resolved {M, N, 0};  // Only alpha electrons to mimic a spin-unresolved case.
    const GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis> linear_expansion_resolved {onv_basis_resolved, linear_expansion_unresolved.coefficients()};

    const auto D_resolved = linear_expansion_resolved.calculate1DM();  // This is the orbital 1-DM, but there are no beta contributions anyways.
    BOOST_CHECK(D_unresolved.isApprox(D_resolved, 1.0e-12));
}

/**
 *  Check some 1-DM values calculated for the SpinUnresolvedONVBasis by comparing them to an equivalent spin-resolved calculation.
 */
BOOST_AUTO_TEST_CASE(calculate1DM_SpinUnresolved_ref_SpinResolved_doubleCheck) {

    // Set up an example linear expansion in a spin-unresolved ONV basis.
    const size_t M = 8;
    const size_t N = 3;
    const GQCP::SpinUnresolvedONVBasis onv_basis {M, N};

    const auto linear_expansion_unresolved = GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis>::Random(onv_basis);
    const auto D_unresolved = linear_expansion_unresolved.calculate1DM();


    // Create an equivalent, spin-resolved linear expansion.
    const GQCP::SpinResolvedONVBasis onv_basis_resolved {M, N, 0};  // Only alpha electrons to mimic a spin-unresolved case.
    const GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis> linear_expansion_resolved {onv_basis_resolved, linear_expansion_unresolved.coefficients()};

    const auto D_resolved = linear_expansion_resolved.calculate1DM();  // This is the orbital 1-DM, but there are no beta contributions anyways.
    BOOST_CHECK(D_unresolved.isApprox(D_resolved, 1.0e-12));
}
