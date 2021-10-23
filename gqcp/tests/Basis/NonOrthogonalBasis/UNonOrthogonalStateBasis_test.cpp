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

#define BOOST_TEST_MODULE "UNonOrthogonalStateBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/NonOrthogonalBasis/UNonOrthogonalStateBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test whether or not the constructor of a non-orthogonal state basis throws when this is expected.
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // This test case is taken from a python prototype from @lelemmen and @johdvos.
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 0.944863, 0);  // H2, 0.5 Angstrom apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize some non-orthogonal "unrestricted states".
    GQCP::SquareMatrix<double> state_1_alpha {4};
    // clang-format off
    state_1_alpha << 0.11206486,  0.0664533 ,  -0.06223693, -0.23572168,
                     0.14935776,  1.61380725,   0.57998371,  0.50570049,
                     0.05498453, -0.00201046,  -0.14404996,  0.42154899,
                     0.06584014, -1.03012307,   0.61860373, -0.38383381;
    // clang-format on
    GQCP::SquareMatrix<double> state_1_beta {4};
    // clang-format off
    state_1_beta <<  0.35722807,   0.01455789, -0.65501355, -0.1817075 ,
                     0.0558135 ,   1.53245558,  0.05524368,  0.65134408,
                     0.26831116,  -0.04871381, -0.47300476,  1.13157038,
                     0.11243095,  -1.29670236,  0.06771168, -0.77975514;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_alpha {4};
    // clang-format off
    state_2_alpha <<  0.26797459,  0.03623556, -0.40680507, -0.33042291,
                      0.19357135,  2.34671211,  0.10421429,  1.4965034 ,
                      0.08067523, -0.00582667, -0.31647872,  0.8164066 ,
                      0.16186438, -2.45695288,  0.45360206, -0.61357143;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_beta {4};
    // clang-format off
    state_2_beta <<  0.340871928  ,  0.0652390350, -0.491818890 , -0.987725549,
                     0.0927480365 ,  2.22620330  ,  0.0325517208,  0.310965988,
                     0.00167619567, -0.0560473064, -0.211679596 ,  0.638032341,
                     0.186865435  , -0.529305142 ,  0.316966132 , -0.376218036;
    // clang-format on
    GQCP::SquareMatrix<double> state_3_alpha {4};
    // clang-format off
    state_3_alpha <<  0.32165359,  0.07409797, -0.61754043, -1.29664545,
                      0.17649753,  0.42795186,  0.14375287,  0.25679393,
                      0.0132411 , -0.00790032, -0.27440902,  0.32891666,
                      0.02341172, -1.21673966,  0.42086602, -0.56974723;
    // clang-format on
    GQCP::SquareMatrix<double> state_3_beta {4};
    // clang-format off
    state_3_beta <<  0.03621239,  0.05062348,  -0.1906441 ,  -1.28517284,
                     0.09895494,  2.0934359 ,   0.20758293,   0.39451284,
                     0.27260356, -0.0417335 ,  -0.48423516,   0.20045802,
                     0.19128949, -0.62531908,   0.19289412,  -1.50848556;
    // clang-format on
    GQCP::SquareMatrix<double> state_4_alpha {3};
    // clang-format off
    state_4_alpha <<  0.04917986,  0.04546257, -0.1548098 ,
                      0.00962112,  0.62442864,  0.46972332,
                      0.35956581, -0.05316816, -0.5599045 ;
    // clang-format on
    GQCP::SquareMatrix<double> state_4_beta {3};
    // clang-format off
    state_4_beta <<  0.28046056,  0.03568668, -0.0270097 , 
                     0.01472352,  0.14086505,  0.58691824,  
                     0.31468179, -0.00655532, -0.1611742 ;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto state1_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_1_alpha}, GQCP::UTransformationComponent<double> {state_1_beta}};
    const auto state2_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_2_alpha}, GQCP::UTransformationComponent<double> {state_2_beta}};
    const auto state3_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_3_alpha}, GQCP::UTransformationComponent<double> {state_3_beta}};
    const auto state4_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_4_alpha}, GQCP::UTransformationComponent<double> {state_4_beta}};

    using NonOrthogonalStateBasisType = GQCP::UNonOrthogonalStateBasis<double>;
    // Check that the constructor doesn't throw an exception when the right dimensions of states are used.
    BOOST_CHECK_NO_THROW(NonOrthogonalStateBasisType basis(std::vector<GQCP::UTransformation<double>> {state1_expansion, state2_expansion, state3_expansion}, S, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()));
    // Check that the constructor throws an exception when the wrong dimensions of states are used.
    BOOST_CHECK_THROW(NonOrthogonalStateBasisType basis_wrong(std::vector<GQCP::UTransformation<double>> {state1_expansion, state2_expansion, state4_expansion}, S, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()), std::invalid_argument);
}

/**
 *  Test whether or not the overlap operator is evaluated correctly over a non-orthogonal state basis.
 */
BOOST_AUTO_TEST_CASE(overlap_and_hamiltonian) {
    // This test case is taken from a python prototype from @lelemmen and @johdvos.
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 0.944863, 0);  // H2, 0.5 Angstrom apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::USpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S = spin_orbital_basis.overlap();

    // Initialize some non-orthogonal "unrestricted states".
    GQCP::SquareMatrix<double> state_1_alpha {4};
    // clang-format off
    state_1_alpha << 0.11206486,  0.0664533 ,  -0.06223693, -0.23572168,
                     0.14935776,  1.61380725,   0.57998371,  0.50570049,
                     0.05498453, -0.00201046,  -0.14404996,  0.42154899,
                     0.06584014, -1.03012307,   0.61860373, -0.38383381;
    // clang-format on
    GQCP::SquareMatrix<double> state_1_beta {4};
    // clang-format off
    state_1_beta <<  0.35722807,   0.01455789, -0.65501355, -0.1817075 ,
                     0.0558135 ,   1.53245558,  0.05524368,  0.65134408,
                     0.26831116,  -0.04871381, -0.47300476,  1.13157038,
                     0.11243095,  -1.29670236,  0.06771168, -0.77975514;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_alpha {4};
    // clang-format off
    state_2_alpha <<  0.26797459,  0.03623556, -0.40680507, -0.33042291,
                      0.19357135,  2.34671211,  0.10421429,  1.4965034 ,
                      0.08067523, -0.00582667, -0.31647872,  0.8164066 ,
                      0.16186438, -2.45695288,  0.45360206, -0.61357143;
    // clang-format on
    GQCP::SquareMatrix<double> state_2_beta {4};
    // clang-format off
    state_2_beta <<  0.340871928  ,  0.0652390350, -0.491818890 , -0.987725549,
                     0.0927480365 ,  2.22620330  ,  0.0325517208,  0.310965988,
                     0.00167619567, -0.0560473064, -0.211679596 ,  0.638032341,
                     0.186865435  , -0.529305142 ,  0.316966132 , -0.376218036;
    // clang-format on
    GQCP::SquareMatrix<double> state_3_alpha {4};
    // clang-format off
    state_3_alpha <<  0.32165359,  0.07409797, -0.61754043, -1.29664545,
                      0.17649753,  0.42795186,  0.14375287,  0.25679393,
                      0.0132411 , -0.00790032, -0.27440902,  0.32891666,
                      0.02341172, -1.21673966,  0.42086602, -0.56974723;
    // clang-format on
    GQCP::SquareMatrix<double> state_3_beta {4};
    // clang-format off
    state_3_beta <<  0.03621239,  0.05062348,  -0.1906441 ,  -1.28517284,
                     0.09895494,  2.0934359 ,   0.20758293,   0.39451284,
                     0.27260356, -0.0417335 ,  -0.48423516,   0.20045802,
                     0.19128949, -0.62531908,   0.19289412,  -1.50848556;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto state1_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_1_alpha}, GQCP::UTransformationComponent<double> {state_1_beta}};
    const auto state2_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_2_alpha}, GQCP::UTransformationComponent<double> {state_2_beta}};
    const auto state3_expansion = GQCP::UTransformation<double> {GQCP::UTransformationComponent<double> {state_3_alpha}, GQCP::UTransformationComponent<double> {state_3_beta}};

    using NonOrthogonalStateBasisType = GQCP::UNonOrthogonalStateBasis<double>;
    NonOrthogonalStateBasisType NOCIbasis {std::vector<GQCP::UTransformation<double>> {state1_expansion, state2_expansion, state3_expansion}, S, molecule.numberOfElectronPairs(), molecule.numberOfElectronPairs()};

    // Evaluate the overlap operator in the NOCI basis.
    const auto overlap_matrix = NOCIbasis.evaluateOverlapOperator();

    // In itialize a reference overlap matrix. Reference data taken from the implementation of @johdvos.
    GQCP::SquareMatrix<double> overlap_reference {3};
    // clang-format off
    overlap_reference <<  0.05573789, 0.07821839, 0.05588782,
                          0.07821839, 0.11900939, 0.07973549,
                          0.05588782, 0.07973549, 0.066858  ;
    // clang-format on

    // Check whether the calculated overlap matches the reference.
    BOOST_CHECK(overlap_matrix.isApprox(overlap_reference, 1e-6));

    // evaluate the Hamiltonian in the NOCI basis.
    const auto sq_ham = spin_orbital_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));
    const auto ham = NOCIbasis.evaluateHamiltonianOperator(sq_ham, GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()));

    // Initialize a reference Hamiltonian. Reference data taken from the implementation of @johdvos.
    GQCP::SquareMatrix<double> hamiltonian_reference {3};
    // clang-format off
    hamiltonian_reference <<  -0.05634827, -0.07880862,  -0.0594536 ,
                              -0.07880862, -0.1082845 ,  -0.08489053,
                              -0.0594536 , -0.08489053,  -0.06239094;
    // clang-format on
    // std::cout << ham << std::endl;
    // Check whether the calculated Hamiltonian matches the reference.
    BOOST_CHECK(ham.isApprox(hamiltonian_reference, 1e-6));
}
