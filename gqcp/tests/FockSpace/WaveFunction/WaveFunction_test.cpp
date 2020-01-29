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
#define BOOST_TEST_MODULE "WaveFunction"

#include <boost/test/unit_test.hpp>

#include "FockSpace/WaveFunction/WaveFunction.hpp"

#include "Basis/transform.hpp"
#include "FockSpace/FockSpace.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"


BOOST_AUTO_TEST_CASE ( shannon_entropy ) {

    // Set up a test Fock space
    GQCP::FockSpace fock_space (8, 3);  // K = 8, N = 3


    // Check the entropy of a Hartree-Fock expansion
    GQCP::WaveFunction hartree_fock_expansion (fock_space, fock_space.HartreeFockExpansion());
    BOOST_CHECK(hartree_fock_expansion.calculateShannonEntropy() < 1.0e-12);  // should be 0


    // Check the maximal entropy (corresponding to a wave function with all equal coefficients different from zero)
    GQCP::WaveFunction constant_expansion (fock_space, fock_space.constantExpansion());
    double reference_entropy = std::log2(fock_space.get_dimension());  // manual derivation
    BOOST_CHECK(std::abs(constant_expansion.calculateShannonEntropy() - reference_entropy) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( transform_wave_function_h3 ) {

    // Produce a wave function, transform it, then pair it against second produced wave function from a transformed basis
    // Create a molecule
    GQCP::Molecule hchain = GQCP::Molecule::HChain(3, 0.742, -1);

    // Create the molecular Hamiltonian for this molecule and basis
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis
    auto K = sq_hamiltonian.dimension();
    auto N_P = hchain.numberOfElectrons()/2;

    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    GQCP::DenseSolverOptions solver_options;

    GQCP::CISolver ci_solver (fci, sq_hamiltonian);
    ci_solver.solve(solver_options);

    // Retrieve the wave function and transform it
    auto wavefunction1 = ci_solver.makeWavefunction();
    GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);

    wavefunction1.basisTransform(U_random);

    // Generate a new wave function by rotating the basis and performing the FCI again
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);
    GQCP::CISolver ci_solver2 (fci, sq_hamiltonian);
    ci_solver2.solve(solver_options);

    auto wavefunction2 = ci_solver2.makeWavefunction();

    // Check if they deviate
    BOOST_CHECK(wavefunction2.isApprox(wavefunction1));
}


BOOST_AUTO_TEST_CASE ( transform_wave_function_h4 ) {

    // Produce a wave function, transform it then pair it, against second produced wave function from a transformed basis
    // Create a molecule
    GQCP::Molecule hchain = GQCP::Molecule::HChain(4, 0.742, 0);

    // Create the molecular Hamiltonian for this molecule and basis
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis
    auto K = sq_hamiltonian.dimension();
    auto N_P = hchain.numberOfElectrons()/2;

    GQCP::ProductFockSpace fock_space (K, N_P, N_P);
    GQCP::FCI fci (fock_space);

    GQCP::DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = 10;
    GQCP::CISolver ci_solver (fci, sq_hamiltonian);
    ci_solver.solve(solver_options);

    // Retrieve the wave function and transform it
    auto wavefunction1 = ci_solver.makeWavefunction();
    GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);

    wavefunction1.basisTransform(U_random);

    // Generate a new wave function by rotating the basis and performing the FCI again.
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);
    GQCP::CISolver ci_solver2 (fci, sq_hamiltonian);
    ci_solver2.solve(solver_options);

    auto wavefunction2 = ci_solver2.makeWavefunction();

    // Check if they deviate
    BOOST_CHECK(wavefunction2.isApprox(wavefunction1));
}


BOOST_AUTO_TEST_CASE ( transform_wave_function_h5 ) {

    // Produce a wave function, transform it then pair it, against second produced wave function from a transformed basis
    // Create a molecule
    GQCP::Molecule hchain = GQCP::Molecule::HChain(5, 0.742, 0);

    // Create the molecular Hamiltonian for this molecule and basis
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (hchain, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, hchain);  // in the Löwdin basis
    auto K = sq_hamiltonian.dimension();
    auto N_B = hchain.numberOfElectrons()/2;
    auto N_A = hchain.numberOfElectrons() - N_B;

    GQCP::ProductFockSpace fock_space (K, N_A, N_B);
    GQCP::FCI fci (fock_space);

    GQCP::DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = 10;

    GQCP::CISolver ci_solver (fci, sq_hamiltonian);
    ci_solver.solve(solver_options);

    // Retrieve the wave function and transform it
    auto wavefunction1 = ci_solver.makeWavefunction();
    GQCP::TransformationMatrix<double> U_random = GQCP::TransformationMatrix<double>::RandomUnitary(K);

    wavefunction1.basisTransform(U_random);

    // Generate a new wave function by rotating the basis and performing the FCI again.
    GQCP::basisRotate(spinor_basis, sq_hamiltonian, U_random);
    GQCP::CISolver ci_solver2 (fci, sq_hamiltonian);
    ci_solver2.solve(solver_options);

    auto wavefunction2 = ci_solver2.makeWavefunction();

    // Check if they deviate
    BOOST_CHECK(wavefunction2.isApprox(wavefunction1));
}


/*
 *  Tests the select method creating a WaveFunctionSelection with the 2 largest coefficients
 */ 
BOOST_AUTO_TEST_CASE ( constructor ) {

    GQCP::ONV onv1 = GQCP::ONV(3, 1, 1);  // 001
    GQCP::ONV onv2 = GQCP::ONV(3, 1, 2);  // 010

    GQCP::VectorX<double> test_coefficients (4);
    test_coefficients << 2.0/7, -5.0/7, 4.0/7, 2.0/7;

    GQCP::SelectedFockSpace selected_fock_space (3, 1, 1);

    selected_fock_space.addConfiguration(onv1.asString(), onv1.asString());  // 001 | 001
    selected_fock_space.addConfiguration(onv1.asString(), onv2.asString());  // 001 | 010
    selected_fock_space.addConfiguration(onv2.asString(), onv1.asString());  // 010 | 001
    selected_fock_space.addConfiguration(onv2.asString(), onv2.asString());  // 010 | 010

    GQCP::WaveFunction wf (selected_fock_space, test_coefficients);

    // Create a WaveFunctionSelection with 2 of the largest coefficients (in this example: configuration 2 and 3)
    GQCP::WaveFunctionSelection selection = wf.select(2);
    const auto& configurations = selection.get_configurations();
    const auto& coefficients = selection.get_coefficients();
    
    // Given the extraction method, the smallest coefficient will come first this is configuration 3: 010 | 001
    BOOST_CHECK(configurations[0].onv_alpha.asString() == "010");
    BOOST_CHECK(configurations[0].onv_beta.asString() == "001");
    BOOST_CHECK(std::abs(coefficients(0) - (4.0/7)) < 1.0e-9);

    // The second largest contribution configuration is number 2: 001 | 010
    BOOST_CHECK(configurations[1].onv_alpha.asString() == "001");
    BOOST_CHECK(configurations[1].onv_beta.asString() == "010");
    BOOST_CHECK(std::abs(coefficients(1) - (-5.0/7)) < 1.0e-9);
}
