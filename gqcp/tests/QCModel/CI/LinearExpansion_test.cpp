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
#define BOOST_TEST_MODULE "LinearExpansion"

#include <boost/test/unit_test.hpp>

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
BOOST_AUTO_TEST_CASE ( reader_test ) {

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
    const auto onv1_alpha = onv1.alphaONV().asString();
    const auto onv1_beta = onv1.betaONV().asString();

    BOOST_CHECK(onv1_alpha == alpha1_ref);
    BOOST_CHECK(onv1_beta == beta1_ref);

    const auto onv2 = linear_expansion.onvBasis().get_configuration(1);
    const auto onv2_alpha = onv2.alphaONV().asString();
    const auto onv2_beta = onv2.betaONV().asString();

    BOOST_CHECK(onv2_alpha == alpha2_ref);
    BOOST_CHECK(onv2_beta == beta2_ref);
}


/**
 *  Check if the calculation of the Shannon entropy is correctly implemented by comparing with a manual calculation.
 */
BOOST_AUTO_TEST_CASE ( shannon_entropy ) {

    // Set up a test spin-resolved ONV basis.
    const GQCP::SpinUnresolvedONVBasis onv_basis (8, 3);  // 8 spinors, 3 electrons


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
BOOST_AUTO_TEST_CASE ( transform_wave_function_h3 ) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(3, 0.742, -1);  // charge -1
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);

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
BOOST_AUTO_TEST_CASE ( transform_wave_function_h4 ) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(4, 0.742);
    const auto N_P = molecule.numberOfElectrons() / 2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_P, N_P);

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
BOOST_AUTO_TEST_CASE ( transform_wave_function_h4 ) {

    // Create the molecular Hamiltonian in the Löwdin basis.
    const auto molecule = GQCP::Molecule::HChain(5, 0.742);
    const auto N_alpha = 3;
    const auto N_beta = 2

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (molecule, "STO-3G");
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Do a dense FCI calculation.
    const GQCP::SpinResolvedONVBasis onv_basis (K, N_alpha, N_beta);

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
