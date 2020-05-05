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

#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Basis/transform.hpp"
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
    const auto onv1 = linear_expansion.onvBasis().get_configuration(0);
    const auto onv1_alpha = onv1.onv(GQCP::Spin::alpha).asString();
    const auto onv1_beta = onv1.onv(GQCP::Spin::beta).asString();

    BOOST_CHECK(onv1_alpha == alpha1_ref);
    BOOST_CHECK(onv1_beta == beta1_ref);

    const auto onv2 = linear_expansion.onvBasis().get_configuration(1);
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
    GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> hartree_fock_expansion {onv_basis, onv_basis.hartreeFockExpansion()};
    BOOST_CHECK(hartree_fock_expansion.calculateShannonEntropy() < 1.0e-12);  // should be 0


    // Check the maximal entropy (corresponding to a wave function with all equal coefficients different from zero)
    GQCP::LinearExpansion<GQCP::SpinUnresolvedONVBasis> constant_expansion {onv_basis, onv_basis.constantExpansion()};
    const double reference_entropy = std::log2(onv_basis.get_dimension());  // manual derivation
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

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation matrix and calculate the transformation of the linear expansion coefficients.
    const GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransformInPlace(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);

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

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_P, N_P};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation matrix and calculate the transformation of the linear expansion coefficients.
    const GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransformInPlace(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);

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

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis {K, N_alpha, N_beta};

    auto environment_direct = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_direct = GQCP::EigenproblemSolver::Dense();

    auto linear_expansion_direct = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_direct, environment_direct).groundStateParameters();


    // Generate a random rotation matrix and calculate the transformation of the linear expansion coefficients.
    const GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);
    linear_expansion_direct.basisTransformInPlace(U_random);


    // Calculate a new linear expansion by rotation the underlying spinor basis and doing another dense calculation, and check if they deviate.
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);

    auto environment_indirect = GQCP::CIEnvironment::Dense(sq_hamiltonian, onv_basis);
    auto solver_indirect = GQCP::EigenproblemSolver::Dense();

    const auto linear_expansion_indirect = GQCP::QCMethod::CI<GQCP::SpinResolvedONVBasis>(onv_basis).optimize(solver_indirect, environment_indirect).groundStateParameters();
    BOOST_CHECK(linear_expansion_direct.isApprox(linear_expansion_indirect, 1.0e-12));
}


// /**
//  *  Check if the projection of |UHF> and |GHF> (the equivalent ONV in a generalized spinor basis) onto |RHF> is the same.
//  */
// BOOST_AUTO_TEST_CASE(overlap_GHF_UHF_on_RHF) {

//     const auto molecule = GQCP::Molecule::HRingFromDistance(4, 1.0);

//     GQCP::RSpinorBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
//     GQCP::TransformationMatrix<double> C {4};  // RHF canonical orbitals for this system (Xeno)
//     // clang-format off
//     C << -0.27745359, -0.8505133,   0.85051937,  2.02075317,
//          -0.27745362, -0.85051937, -0.8505133,  -2.02075317,
//          -0.27745359,  0.8505133,  -0.85051937,  2.02075317,
//          -0.27745362,  0.85051937,  0.8505133,  -2.02075317;
//     // clang-format on
//     r_spinor_basis.transform(C);


//     // GQCP::USpinorBasis<double, GQCP::GTOShell> u_spinor_basis {molecule, "STO-3G"};
//     GQCP::TransformationMatrix<double> C_alpha {4};  // UHF alpha canonical orbitals for this system (Xeno), triplet
//     // clang-format off
//     C_alpha << -1.75646828e-01, -1.20606646e-06,  1.20281173e+00,  2.03213486e+00,
//                -3.78560533e-01, -1.20281173e+00, -1.20606647e-06, -2.00427438e+00,
//                -1.75646828e-01,  1.20606646e-06, -1.20281173e+00,  2.03213486e+00,
//                -3.78560533e-01,  1.20281173e+00,  1.20606646e-06, -2.00427438e+00;
//     // clang-format on

//     GQCP::TransformationMatrix<double> C_beta {4};  // UHF alpha canonical orbitals for this system (Xeno), triplet
//     // clang-format off
//     C_beta << -3.78560533e-01,  1.20281173e+00,  1.21724557e-06,  2.00427438e+00,
//               -1.75646828e-01,  1.21724558e-06, -1.20281173e+00, -2.03213486e+00,
//               -3.78560533e-01, -1.20281173e+00, -1.21724558e-06,  2.00427438e+00,
//               -1.75646828e-01, -1.21724558e-06,  1.20281173e+00, -2.03213486e+00;
//     // clang-format on
//     // u_spinor_basis.transform(C_alpha, C_beta);


//     const auto u_spinor_basis = GQCP::USpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);


//     const auto uhf_onv = GQCP::SpinResolvedONV::UHF(4, 2, 2);
//     const auto linear_expansion = GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>::FromONVProjection(uhf_onv, r_spinor_basis, u_spinor_basis);

//     std::cout << linear_expansion.coefficients() << std::endl;
// }
