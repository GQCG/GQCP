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

#define BOOST_TEST_MODULE "NOCI_test"

#include <boost/test/unit_test.hpp>

#include "Basis/BiorthogonalBasis/ULowdinPairingBasis.hpp"
#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Mathematical/Optimization/Eigenproblem/GeneralizedEigenproblemSolver.hpp"
#include "Molecule/Molecule.hpp"
#include "QCMethod/NOCI/NOCI.hpp"
#include "QCMethod/NOCI/NOCIEnvironment.hpp"
#include "QCModel/NOCI/NOCIExpansion.hpp"
#include "Utilities/complex.hpp"

/**
 *  Test the complex unrestricted NOCI model 1DM calculation.
 */
BOOST_AUTO_TEST_CASE(NOCI_unrestricted_crash_test) {
    // This test case is taken from a python prototype from @lelemmen and @johdvos.
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 3.023561581760001);  // H2, 2.5 bohr apart.

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "restricted states".
    GQCP::SquareMatrix<double> state_1_a {3};
    GQCP::SquareMatrix<double> state_1_b {3};
    // clang-format off
    state_1_a << 3.01612395e-01, -7.27196084e-15,  1.18334666e+00,
               4.60055205e-01, -9.96503166e-01, -5.35359699e-01,
               4.60055205e-01,  9.96503166e-01, -5.35359699e-01;

    state_1_b << 6.38858962e-01, -1.04073944e+00, -2.27595721e-15,
                 2.80666711e-01,  6.47678158e-01, -9.96503166e-01,
                 2.80666711e-01,  6.47678158e-01,  9.96503166e-01;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_a {3};
    GQCP::SquareMatrix<double> state_2_b {3};
    // clang-format off
    state_2_a <<  0.30161239,  0.63885896,  1.18334666,
                  0.46005521,  0.28066671, -0.5353597,
                  0.46005521,  0.28066671, -0.5353597;

    state_2_b << -7.27196084e-15, -1.04073944e+00, -2.27595721e-15,
                 -9.96503166e-01,  6.47678158e-01, -9.96503166e-01,
                  9.96503166e-01,  6.47678158e-01,  9.96503166e-01;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_1_a}, GQCP::UTransformationComponent<double> {state_1_b}};
    const auto basis_state_2 = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_2_a}, GQCP::UTransformationComponent<double> {state_2_b}};

    const auto lowdin_pairing_basis = GQCP::ULowdinPairingBasis<double>(basis_state_1, basis_state_2, S, 2, 1);
    // std::cout << lowdin_pairing_basis.numberOfZeroOverlaps() << std::endl;
    // std::cout << lowdin_pairing_basis.zeroOverlapIndices()[0].second << std::endl;
    // std::cout << "===============" << std::endl;

    // Create a vector out of these three basis states.
    std::vector<GQCP::UTransformation<double>> basis_vector {basis_state_1, basis_state_2};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::UNonOrthogonalStateBasis<double> {basis_vector, S, 2, 1};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto sq_hamiltonian = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Create a dense solver and corresponding environment and put them together in the QCMethod.
    auto environment = GQCP::NOCIEnvironment::Dense(sq_hamiltonian, NOS_basis, molecule);
    auto solver = GQCP::GeneralizedEigenproblemSolver::Dense<double>();

    // We do not specify the number of states, meaning we only request the ground state.
    const auto NOCI_model = GQCP::QCMethod::NOCI<double, GQCP::UNonOrthogonalStateBasis<double>>(NOS_basis).optimize(solver, environment);

    // In this test case, the alpha and beta 1DM's should be equal. We check this fact, to confirm that all aspects of the complex NOCI method in an unrestricted basis work as expected.
    std::cout << NOCI_model.groundStateEnergy();
    // BOOST_CHECK(std::abs(-1.35 - NOCI_model.groundStateEnergy()) < 1e-6);
}
