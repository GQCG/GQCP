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

#define BOOST_TEST_MODULE "RNonOrthogonalStateBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/NonOrthogonalBasis/RNonOrthogonalStateBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test whether or not the constructor of a non-orthogonal state basis throws when this is expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize some non-orthogonal "Restricted states".
    GQCP::SquareMatrix<double> state_1 {4};
    // clang-format off
    state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,
               -0.07874922,  0.15086478, -0.68085546,  0.77423311,
               -0.24580188,  0.26338108,  0.09556297, -0.12178159,
               -0.38944259,  0.4101685 ,  0.45214166, -0.58335985;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {4};
    // clang-format off
    state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,
                0.36593356, -0.28669343, -0.84796858, -0.13503625,
                0.25853403,  0.14539669,  0.17176599, -0.01126146,
                0.36597032,  0.28670189,  0.847938  , -0.13501526;
    // clang-format on
    GQCP::SquareMatrix<double> state_3 {3};
    // clang-format off
    state_3 <<  0.25851329, -0.14539151, -0.17177142,
                0.36593356, -0.28669343, -0.84796858,
                0.25853403,  0.14539669,  0.17176599;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1_expansion = GQCP::RTransformation<double> {state_1};
    const auto state2_expansion = GQCP::RTransformation<double> {state_2};
    const auto state3_expansion = GQCP::RTransformation<double> {state_3};

    using NonOrthogonalStateBasisType = GQCP::RNonOrthogonalStateBasis<double>;
    // Check that the constructor doesn't throw an exception when the right dimensions of states are used.
    BOOST_CHECK_NO_THROW(NonOrthogonalStateBasisType basis(std::vector<GQCP::RTransformation<double>> {state1_expansion, state2_expansion}, S, molecule.numberOfElectronPairs()));
    // Check that the constructor throws an exception when the wrong dimensions of states are used.
    BOOST_CHECK_THROW(NonOrthogonalStateBasisType basis_wrong(std::vector<GQCP::RTransformation<double>> {state1_expansion, state3_expansion}, S, molecule.numberOfElectronPairs()), std::invalid_argument);
}


/**
 *  Test whether the non-orthogonal state basis evaluates matrix elements correctly.
 */
BOOST_AUTO_TEST_CASE(operator_evaluations) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "restricted states".
    GQCP::SquareMatrix<double> state_1 {4};
    // clang-format off
    state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,
               -0.07874922,  0.15086478, -0.68085546,  0.77423311,
               -0.24580188,  0.26338108,  0.09556297, -0.12178159,
               -0.38944259,  0.4101685 ,  0.45214166, -0.58335985;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {4};
    // clang-format off
    state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,
                0.36593356, -0.28669343, -0.84796858, -0.13503625,
                0.25853403,  0.14539669,  0.17176599, -0.01126146,
                0.36597032,  0.28670189,  0.847938  , -0.13501526;
    // clang-format on
    GQCP::SquareMatrix<double> state_3 {4};
    // clang-format off
    state_3 <<  -0.265842  ,  0.17716735, -0.15969328, -0.00308706,
                -0.36278694,  0.36406651, -0.80340861, -0.13144475,
                -0.26584976, -0.17716927,  0.15969112, -0.00308558,
                -0.36280035, -0.36406982,  0.80339638, -0.13143372;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::RTransformation<double> {state_1};
    const auto basis_state_2 = GQCP::RTransformation<double> {state_2};
    const auto basis_state_3 = GQCP::RTransformation<double> {state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::RTransformation<double>> basis_vector {basis_state_1, basis_state_2, basis_state_3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::RNonOrthogonalStateBasis<double> {basis_vector, S, molecule.numberOfElectronPairs()};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto hamiltonian_AO = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Initialize the reference overlap matrix.
    GQCP::SquareMatrix<double> reference_overlap {3};
    // clang-format off
    reference_overlap <<   0.181724, 0.350152,  0.35311,
                           0.350152, 0.874972, 0.882449,
                           0.35311 , 0.882449, 0.890075;
    // clang-format on

    // Initialize the reference core hamiltonian.
    GQCP::SquareMatrix<double> reference_core_hamiltonian {3};
    // clang-format off
    reference_core_hamiltonian << -0.339766,  -0.67903,  -0.685209,
                                  -0.67903 ,  -1.69767,  -1.7132  ,
                                  -0.685209,   -1.7132,  -1.72894 ;
    // clang-format on

    // Initialize the reference two electron contribution.
    GQCP::SquareMatrix<double> reference_two_electron {3};
    // clang-format off
    reference_two_electron << 0.103774, 0.185861, 0.188009,
                              0.185861, 0.446029, 0.451139,
                              0.188009, 0.451139, 0.456353;
    // clang-format on
    
    BOOST_CHECK(NOS_basis.evaluateOverlapOperator().isApprox(reference_overlap, 1e-6));
    BOOST_CHECK(NOS_basis.evaluateOneElectronOperator(hamiltonian_AO.core()).isApprox(reference_core_hamiltonian, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateTwoElectronOperator(hamiltonian_AO.twoElectron()).isApprox(reference_two_electron, 1e-5));
}
