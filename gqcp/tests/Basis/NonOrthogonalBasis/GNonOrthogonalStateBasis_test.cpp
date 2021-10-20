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

#define BOOST_TEST_MODULE "GNonOrthogonalStateBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test whether or not the constructor of a non-orthogonal state basis throws when this is expected.
 *  This test checks the `GNonOrthogonalStateBasis` case, however as both the `R` and the `G` case test the functionality of the `SimpleNonOrthogonalStateBasis`, the restricted case won't be tested separately.
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // This test case is taken from a python prototype from @lelemmen and @johdvos.
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize two non-orthogonal "Generalised states".
    GQCP::SquareMatrix<double> C_bra {8};
    // clang-format off
    C_bra << -0.07443693,  0.12036042, -0.13557067,  0.15517005,  0.13315100,  -0.03074946, -0.92997876, -0.93718779,
             -0.07874922,  0.15086478, -0.68085546,  0.77423311,  0.08712373,   0.25266303,  0.848079  ,  0.89108911,
             -0.24580188,  0.26338108,  0.09556297, -0.12178159,  0.91319299,   0.90475733, -0.03994767,  0.12839983,
             -0.38944259,  0.4101685 ,  0.45214166, -0.58335985, -0.90125958,  -0.93270816, -0.16410814, -0.32074956,
             -0.26338108, -0.24580188, -0.12178159, -0.09556297, -0.90475733,   0.91319299, -0.12839983, -0.03994767,
             -0.4101685 , -0.38944259, -0.58335985, -0.45214166,  0.93270817,  -0.90125958,  0.32074956, -0.16410814,
             -0.12036042, -0.07443693,  0.15517005,  0.13557067,  0.03074946,   0.13315101,  0.93718779, -0.92997876,
             -0.15086478, -0.07874922,  0.77423311,  0.68085546, -0.25266303,   0.08712373, -0.89108911,  0.84807900;
    // clang-format on
    GQCP::SquareMatrix<double> C_ket_right {8};
    // clang-format off
    C_ket_right <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
                    0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
                    0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
                    0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
                    0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
                    0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
                   -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
                   -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    GQCP::SquareMatrix<double> C_ket_wrong {7};
    // clang-format off
    C_ket_wrong <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389,
                    0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,
                    0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  ,
                    0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,
                    0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,
                    0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639,
                   -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<double> {C_bra};
    const auto ket_expansion_right = GQCP::GTransformation<double> {C_ket_right};
    const auto ket_expansion_wrong = GQCP::GTransformation<double> {C_ket_wrong};

    using NonOrthogonalStateBasisType = GQCP::GNonOrthogonalStateBasis<double>;
    // Check that the constructor doesn't throw an exception when the right dimensions of states are used.
    BOOST_CHECK_NO_THROW(NonOrthogonalStateBasisType basis(std::vector<GQCP::GTransformation<double>> {bra_expansion, ket_expansion_right}, S, molecule.numberOfElectrons()));
    // Check that the constructor throws an exception when the wrong dimensions of states are used.
    BOOST_CHECK_THROW(NonOrthogonalStateBasisType basis_wrong(std::vector<GQCP::GTransformation<double>> {bra_expansion, ket_expansion_wrong}, S, molecule.numberOfElectrons()), std::invalid_argument);
}


/**
 *  Test whether the non-orthogonal state basis evaluates the Hamiltonian correctly. All other evaluations (one- and two-electron operators) are used within the Hamiltonian functionality.
 *  This test checks the `GNonOrthogonalStateBasis` case, however as both the `R` and the `G` case test the functionality of the `SimpleNonOrthogonalStateBasis`, the restricted case won't be tested separately.
 */
BOOST_AUTO_TEST_CASE(operator_evaluations) {
    // This test case is taken from a python prototype from @lelemmen and @johdvos.
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize two non-orthogonal "Generalised states".
    GQCP::SquareMatrix<double> basis_state_1 {8};
    // clang-format off
    basis_state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,  0.13315100,  -0.03074946, -0.92997876, -0.93718779,
                     -0.07874922,  0.15086478, -0.68085546,  0.77423311,  0.08712373,   0.25266303,  0.848079  ,  0.89108911,
                     -0.24580188,  0.26338108,  0.09556297, -0.12178159,  0.91319299,   0.90475733, -0.03994767,  0.12839983,
                     -0.38944259,  0.4101685 ,  0.45214166, -0.58335985, -0.90125958,  -0.93270816, -0.16410814, -0.32074956,
                     -0.26338108, -0.24580188, -0.12178159, -0.09556297, -0.90475733,   0.91319299, -0.12839983, -0.03994767,
                     -0.4101685 , -0.38944259, -0.58335985, -0.45214166,  0.93270817,  -0.90125958,  0.32074956, -0.16410814,
                     -0.12036042, -0.07443693,  0.15517005,  0.13557067,  0.03074946,   0.13315101,  0.93718779, -0.92997876,
                     -0.15086478, -0.07874922,  0.77423311,  0.68085546, -0.25266303,   0.08712373, -0.89108911,  0.84807900;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_2 {8};
    // clang-format off
    basis_state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
                      0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
                      0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
                      0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
                      0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
                      0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
                     -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
                     -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_3 {8};
    // clang-format off
    basis_state_3 <<  -0.265842  ,  0.17716735, -0.15969328, -0.00308706,  0.84741422, -0.81942223,  0.41366608, -0.36889812,
                      -0.36278694,  0.36406651, -0.80340861, -0.13144475, -0.65133121,  1.03324716, -0.44428872,  0.24605534,
                      -0.26584976, -0.17716927,  0.15969112, -0.00308558,  0.84933355,  0.81745234, -0.41355441, -0.36897416,
                      -0.36280035, -0.36406982,  0.80339638, -0.13143372, -0.65375477, -1.0317334 ,  0.4442138 ,  0.24613641,
                      -0.09736842,  0.22594293,  0.06676532,  0.17043252,  0.3439281 , -0.39253422, -0.84701679,  0.86052195,
                      -0.15038318,  0.32968947,  0.23916148,  0.91374992, -0.47007004,  0.32909412,  0.60131983, -0.98100354,
                       0.0973641 ,  0.22593416,  0.06676626, -0.17043417, -0.34482098, -0.39170837, -0.84682078, -0.86073623,
                       0.15037797,  0.32967347,  0.2391655 , -0.91376119,  0.47081908,  0.32796701,  0.60109455,  0.98115454;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<double> {basis_state_1};
    const auto ket_expansion_right = GQCP::GTransformation<double> {basis_state_2};
    const auto ket_expansion_wrong = GQCP::GTransformation<double> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {basis_state_1, basis_state_2, basis_state_3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, molecule.numberOfElectrons()};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto hamiltonian_AO = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // We can now evaluate this Hamiltonian operator in the non-orthogonal basis.
    // The Hamiltonian evaluated in non-orthogonal basis requires the nuclear repuslion to be taken into account.
    const auto non_orthogonal_hamiltonian = NOS_basis.evaluateHamiltonianOperator(hamiltonian_AO, GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()));

    // Set up a reference Hamiltonian.
    // The dimension of the Hamiltonian is equal to the number of basis states used.
    GQCP::SquareMatrix<double> reference_hamiltonian {3};
    // clang-format off
    reference_hamiltonian <<  -1.03611672, -0.83337779, -0.78876307,
                              -0.83337779, -1.03337473, -1.02413064,
                              -0.78876307, -1.02413064, -1.0270364 ;
    // clang-format on

    // Compare the reference with the evaluated Hamiltonian.
    BOOST_CHECK(non_orthogonal_hamiltonian.isApprox(reference_hamiltonian, 1e-6));
}
