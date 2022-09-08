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
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "6-31G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "Generalised states".
    GQCP::SquareMatrix<double> state_1 {8};
    // clang-format off
    state_1 << -0.07443693,  0.12036042, -0.13557067,  0.15517005,  0.13315100,  -0.03074946, -0.92997876, -0.93718779,
               -0.07874922,  0.15086478, -0.68085546,  0.77423311,  0.08712373,   0.25266303,  0.848079  ,  0.89108911,
               -0.24580188,  0.26338108,  0.09556297, -0.12178159,  0.91319299,   0.90475733, -0.03994767,  0.12839983,
               -0.38944259,  0.4101685 ,  0.45214166, -0.58335985, -0.90125958,  -0.93270816, -0.16410814, -0.32074956,
               -0.26338108, -0.24580188, -0.12178159, -0.09556297, -0.90475733,   0.91319299, -0.12839983, -0.03994767,
               -0.4101685 , -0.38944259, -0.58335985, -0.45214166,  0.93270817,  -0.90125958,  0.32074956, -0.16410814,
               -0.12036042, -0.07443693,  0.15517005,  0.13557067,  0.03074946,   0.13315101,  0.93718779, -0.92997876,
               -0.15086478, -0.07874922,  0.77423311,  0.68085546, -0.25266303,   0.08712373, -0.89108911,  0.84807900;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {8};
    // clang-format off
    state_2 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
                0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
                0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
                0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
                0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
                0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
               -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
               -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    GQCP::SquareMatrix<double> state_3 {7};
    // clang-format off
    state_3 <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389,
                0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,
                0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  ,
                0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,
                0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,
                0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639,
               -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1_expansion = GQCP::GTransformation<double> {state_1};
    const auto state2_expansion = GQCP::GTransformation<double> {state_2};
    const auto state3_expansion = GQCP::GTransformation<double> {state_3};

    using NonOrthogonalStateBasisType = GQCP::GNonOrthogonalStateBasis<double>;
    // Check that the constructor doesn't throw an exception when the right dimensions of states are used.
    BOOST_CHECK_NO_THROW(NonOrthogonalStateBasisType basis(std::vector<GQCP::GTransformation<double>> {state1_expansion, state2_expansion}, S, molecule.numberOfElectrons()));
    // Check that the constructor throws an exception when the wrong dimensions of states are used.
    BOOST_CHECK_THROW(NonOrthogonalStateBasisType basis_wrong(std::vector<GQCP::GTransformation<double>> {state1_expansion, state3_expansion}, S, molecule.numberOfElectrons()), std::invalid_argument);
}


/**
 *  Test whether the non-orthogonal state basis evaluates the Hamiltonian correctly. All other evaluations (one- and two-electron operators) are used within the Hamiltonian functionality.
 */
BOOST_AUTO_TEST_CASE(operator_evaluations_3_states) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.70075338974);  // H3, 0.9Å apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<double> basis_state_1 {6};
    // clang-format off
    basis_state_1 << 0.5666496, -1.0626127595,   0.0        ,  -0.4614557058, 0.0          ,   0.0         ,
                     0.5666496,  1.0626128447,   0.0        ,  -0.4614554991, 0.0          ,   0.0         ,
                     0.0      , -1.254e-07   ,   0.0        ,  1.28964296040, 0.0          ,   0.0         ,
                     0.0      ,  0.0         , -7.1e-09     ,  0.0          , -0.7307757084, -1.0626127372 ,
                     0.0      ,  0.0         ,  7.1e-09     ,  0.0          , -0.7307755196,  1.0626128671 ,
                     0.0      ,  0.0         , -1.0000000046,  0.0          ,  0.8143580024, -6.57e-080    ;
    // clang-format on
    GQCP::SquareMatrix<double> basis_state_2 {6};
    // clang-format off
    basis_state_2 <<  0.5666495639, -1.0626128318,   0.0         ,  0.46145553040000004,   0.0         ,  0.0         ,
                      0.0         ,  8.74e-08    ,   0.0         , -1.2896429604000001 ,   0.0         ,  0.0         ,
                      0.566649558 ,  1.0626127724,   0.0         ,  0.4614556745       ,   0.0         ,  0.0         ,
                      0.0         ,  0.0         ,  -4.9e-09     ,  0.0                ,   0.7307755482,  1.0626128474,
                      0.0         ,  0.0         ,   1.0000000046,  0.0                ,  -0.8143580024, -4.58e-08    ,
                      0.0         ,  0.0         ,   4.9e-09     ,  0.0                ,   0.7307756798, -1.0626127568;

    // clang-format on
    GQCP::SquareMatrix<double> basis_state_3 {6};
    // clang-format off
    basis_state_3 <<  0.0         ,  1.071e-07   ,   0.0                  ,  -1.2896429604,   0.0         ,  0.0         ,
                      0.5666495646, -1.0626128385,   0.0                  ,   0.4614555142,   0.0         ,  0.0         ,
                      0.5666495574,  1.0626127657,   0.0                  ,   0.4614556907,   0.0         ,  0.0         ,
                      0.0         ,  0.0         ,  -1.0000000046         ,   0.0         ,  -0.8143580024, -5.61e-08    ,
                      0.0         ,  0.0         ,   6.000000000000001e-09,   0.0         ,   0.7307755334,  1.0626128576,
                      0.0         ,  0.0         ,  -6.000000000000001e-09,   0.0         ,   0.7307756946, -1.0626127467;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1 = GQCP::GTransformation<double> {basis_state_1};
    const auto state2 = GQCP::GTransformation<double> {basis_state_2};
    const auto state3 = GQCP::GTransformation<double> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {state1, state2, state3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, molecule.numberOfElectrons()};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto hamiltonian_AO = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // We can now evaluate this Hamiltonian operator in the non-orthogonal basis.
    // The Hamiltonian evaluated in non-orthogonal basis requires the nuclear repuslion to be taken into account.
    const auto non_orthogonal_hamiltonian = NOS_basis.evaluateHamiltonianOperator(hamiltonian_AO, GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()));

    // Initialize the reference overlap matrix.
    GQCP::SquareMatrix<double> reference_overlap {3};
    // clang-format off
    reference_overlap <<   1.0     ,  -0.199371, -0.199371,
                          -0.199371,   1.0     , -0.199371,
                          -0.199371,  -0.199371,  1.0     ;
    // clang-format on

    // Initialize the reference core hamiltonian.
    GQCP::SquareMatrix<double> reference_core_hamiltonian {3};
    // clang-format off
    reference_core_hamiltonian << -4.39452,  1.0116 ,   1.0116 ,
                                   1.0116 , -4.39452,   1.0116 ,
                                   1.0116 ,   1.0116,  -4.39452;
    // clang-format on

    // Initialize the reference two electron contribution.
    GQCP::SquareMatrix<double> reference_two_electron {3};
    // clang-format off
    reference_two_electron << 1.45235,   -0.322416,   -0.322416,
                             -0.322416,   1.45235 ,   -0.322416,
                             -0.322416,  -0.322416,    1.45235 ;
    // clang-format on

    BOOST_CHECK(NOS_basis.evaluateOverlapOperator().isApprox(reference_overlap, 1e-6));
    BOOST_CHECK(NOS_basis.evaluateOneElectronOperator(hamiltonian_AO.core()).isApprox(reference_core_hamiltonian, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateTwoElectronOperator(hamiltonian_AO.twoElectron()).isApprox(reference_two_electron, 1e-5));
}


/**
 *  Test whether the non-orthogonal state basis evaluates the Hamiltonian correctly. All other evaluations (one- and two-electron operators) are used within the Hamiltonian functionality.
 */
BOOST_AUTO_TEST_CASE(operator_evaluations_complex) {
    // This is a self implemented test case, with rotated H3 states in STO-3G.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.88973, 0);  // H3, 1 Angstrom apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<GQCP::complex> basis_state_1 {6};
    using namespace std::complex_literals;
    // clang-format off
    basis_state_1 << -0.40144765 + 0.0i       ,  0.55646711 - 0.0i       ,  0.0        + 0.0i       ,  -0.72745044 + 0.0i      ,   0.0        + 0.0i       ,  0.0        + 0.0i      ,
                      0.02868493 + 0.71242345i, -0.31531155 - 0.19005498i,  0.0        + 0.0i       ,  -0.2570292  - 0.5385385i,   0.0        + 0.0i       ,  0.0        + 0.0i      ,
                      0.50792928 - 0.26921669i,  0.38680437 - 0.63654102i,  0.0        + 0.0i       ,   0.01558441 - 0.3383567i,   0.0        + 0.0i       ,  0.0        + 0.0i      ,
                      0.0        + 0.0i       ,  0.0        + 0.0i        , -0.40144765 + 0.0i       ,   0.0        + 0.0i       ,   0.55646711 - 0.0i       , -0.72745044 + 0.0i      ,
                      0.0        + 0.0i       ,  0.0        + 0.0i        ,  0.02868493 + 0.71242345i,   0.0        + 0.0i       ,  -0.31531155 - 0.19005498i, -0.2570292  - 0.5385385i,
                      0.0        + 0.0i       ,  0.0        + 0.0i        ,  0.50792928 - 0.26921669i,   0.0        + 0.0i       ,   0.38680437 - 0.63654102i,  0.01558441 - 0.3383567i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> basis_state_2 {6};
    // clang-format off
    basis_state_2 <<  -0.61461988 + 0.0i       ,  0.56165315 + 0.0i       ,  0.0        + 0.0i       , -0.5538846  + 0.0i       ,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                      -0.06381077 + 0.37763728i, -0.62473282 + 0.38125827i,  0.0        + 0.0i       , -0.56268722 - 0.03244082i,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                       0.68058226 + 0.1112136i ,  0.2804094  - 0.2650799i ,  0.0        + 0.0i       , -0.47086806 - 0.39220634i,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                       0.0        + 0.0i       ,  0.0        + 0.0i       , -0.61461988 + 0.0i       ,  0.0        + 0.0i       ,  0.56165315 + 0.0i       , -0.5538846  + 0.0i       ,
                       0.0        + 0.0i       ,  0.0        + 0.0i       , -0.06381077 + 0.37763728i,  0.0        + 0.0i       , -0.62473282 + 0.38125827i, -0.56268722 - 0.03244082i,
                       0.0        + 0.0i       ,  0.0        + 0.0i       ,  0.68058226 + 0.1112136i ,  0.0        + 0.0i       ,  0.2804094  - 0.2650799i , -0.47086806 - 0.39220634i;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> basis_state_3 {6};
    // clang-format off
    basis_state_3 <<  -0.76945889 + 0.0i       , -0.43823384 + 0.0i       ,  0.0        + 0.0i       , -0.46463332 + 0.0i       ,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                      -0.12309236 + 0.29558626i,  0.69338911 - 0.02945486i,  0.0        + 0.0i       , -0.45014435 - 0.46172616i,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                       0.37337164 + 0.40743548i, -0.56583894 - 0.07823902i,  0.0        + 0.0i       , -0.08463526 - 0.6009424i ,  0.0        + 0.0i       ,  0.0        + 0.0i       ,
                       0.0        + 0.0i       ,  0.0        + 0.0i       , -0.76945889 + 0.0i       ,  0.0        + 0.0i       , -0.43823384 + 0.0i       , -0.46463332 + 0.0i       ,
                       0.0        + 0.0i       ,  0.0        + 0.0i       , -0.12309236 + 0.29558626i,  0.0        + 0.0i       ,  0.69338911 - 0.02945486i, -0.45014435 - 0.46172616i,
                       0.0        + 0.0i       ,  0.0        + 0.0i       ,  0.37337164 + 0.40743548i,  0.0        + 0.0i       , -0.56583894 - 0.07823902i, -0.08463526 - 0.6009424i ;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto state1 = GQCP::GTransformation<GQCP::complex> {basis_state_1};
    const auto state2 = GQCP::GTransformation<GQCP::complex> {basis_state_2};
    const auto state3 = GQCP::GTransformation<GQCP::complex> {basis_state_3};

    // Create a vector out of these three basis states.
    std::vector<GQCP::GTransformation<GQCP::complex>> basis_vector {state1, state2, state3};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<GQCP::complex> {basis_vector, S, molecule.numberOfElectrons()};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto hamiltonian_AO = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // We can now evaluate this Hamiltonian operator in the non-orthogonal basis.
    // The Hamiltonian evaluated in non-orthogonal basis requires the nuclear repuslion to be taken into account.
    const auto non_orthogonal_hamiltonian = NOS_basis.evaluateHamiltonianOperator(hamiltonian_AO, GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()));

    /// Initialize the reference overlap matrix.
    GQCP::SquareMatrix<GQCP::complex> reference_overlap {3};
    // clang-format off
    reference_overlap <<   0.351783  + 5.61784e-17i,  0.138083  + 0.0987121i  ,  0.0976799 - 0.246505i   ,
                           0.138083  - 0.0987121i  ,  0.204352  - 4.00627e-17i, -0.0959736 - 0.211324i   ,
                           0.0976799 + 0.246505i   , -0.0959736 + 0.211324i   ,  0.4162    - 9.83463e-17i;
    // clang-format on

    // Initialize the reference core hamiltonian.
    GQCP::SquareMatrix<GQCP::complex> reference_core_hamiltonian {3};
    // clang-format off
    reference_core_hamiltonian << -1.34571  - 1.74712e-16i,  -0.492403 - 0.375284i  ,  -0.429008 + 0.956919i   ,
                                  -0.492403 + 0.375284i   ,  -0.725901 + 1.2722e-16i,   0.343079 + 0.793144i   ,
                                  -0.429008 - 0.956919i   ,   0.343079 - 0.793144i  ,  -1.64854  + 2.02519e-16i;
    // clang-format on

    // Initialize the reference two electron contribution.
    GQCP::SquareMatrix<GQCP::complex> reference_two_electron {3};
    // clang-format off
    reference_two_electron << 0.561229 + 4.36825e-17i,   0.22163  + 0.152087i   ,   0.157533 - 0.386507i   ,
                              0.22163  - 0.152087i   ,   0.32012  - 3.70326e-17i,  -0.156274 - 0.324873i   ,
                              0.157533 + 0.386507i   ,  -0.156274 + 0.324873i   ,   0.658129 - 4.34946e-18i;
    // clang-format on

    BOOST_CHECK(NOS_basis.evaluateOverlapOperator().isApprox(reference_overlap, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateOneElectronOperator(hamiltonian_AO.core()).isApprox(reference_core_hamiltonian, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateTwoElectronOperator(hamiltonian_AO.twoElectron()).isApprox(reference_two_electron, 1e-5));
}


/**
 *  Test the GNOCI in the case where one zero overlap value is present. Test case should yield the same results as the UNonOrthogonalStateBasis.
 */
BOOST_AUTO_TEST_CASE(NOCI_from_unrestricted_one_zero) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    // It was for H3, at 1.6Å internuclear distance for the STO-3G basis set.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 3.023561581760001);  // H3, 1.6Å apart.

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap();

    // Initialize two non-orthogonal "generalized states". These are constructed from the unrestricted case.
    // Columns 2 and 3 are swapped to make sure the correct occupied orbitals are selected by the algorithm.
    GQCP::SquareMatrix<double> state_1 {6};
    // clang-format off
    state_1 << 0.639638709, -0.801711730,  0.0,  0.189534772,  0.0        ,  0.0        ,
               0.0        ,  0.0        ,  0.0, -1.04297785 ,  0.0        ,  0.0        , 
               0.639638715,  0.801711722,  0.0,  0.189534786,  0.0        ,  0.0        ,
               0.0        ,  0.0        ,  0.0,  0.0        ,  0.667128998, -0.801711735,
               0.0        ,  0.0        ,  1.0,  0.0        , -0.296315367,  0.0        ,
               0.0        ,  0.0        ,  0.0,  0.0        ,  0.667129020,  0.801711717;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {6};
    // clang-format off
    state_2 <<  0.0,  0.667128998,  0.0        , -0.801711735,  0.0        ,  0.0        ,
                1.0, -0.296315367,  0.0        ,  0.0        ,  0.0        ,  0.0        ,
                0.0,  0.667129020,  0.0        ,  0.801711717,  0.0        ,  0.0        ,
                0.0,  0.0        ,  0.639638709,  0.0        , -0.801711730,  0.189534772,
                0.0,  0.0        ,  0.0        ,  0.0        ,  0.0        , -1.04297785 ,
                0.0,  0.0        ,  0.639638715,  0.0        ,  0.801711722,  0.189534786;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::GTransformation<double> {state_1};
    const auto basis_state_2 = GQCP::GTransformation<double> {state_2};

    // Create a vector out of these two basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {basis_state_1, basis_state_2};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, 3};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Initialize the reference overlap matrix.
    GQCP::SquareMatrix<double> reference_overlap {2};
    // clang-format off
    reference_overlap <<  1.0,  0.0,
                          0.0,  1.0;
    // clang-format on

    // Initialize the reference core hamiltonian.
    GQCP::SquareMatrix<double> reference_core_hamiltonian {2};
    // clang-format off
    reference_core_hamiltonian << -3.32916    ,  -1.53075e-09,
                                  -1.53075e-09,  -3.47686    ;
    // clang-format on

    // Initialize the reference two electron contribution.
    GQCP::SquareMatrix<double> reference_two_electron {2};
    // clang-format off
    reference_two_electron << 0.96323    ,  5.61804e-10,
                              5.61804e-10,  1.18816    ;
    // clang-format on

    BOOST_CHECK(NOS_basis.evaluateOverlapOperator().isApprox(reference_overlap, 1e-6));
    BOOST_CHECK(NOS_basis.evaluateOneElectronOperator(sq_hamiltonian.core()).isApprox(reference_core_hamiltonian, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateTwoElectronOperator(sq_hamiltonian.twoElectron()).isApprox(reference_two_electron, 1e-5));
}


/**
 *  Test the GNOCI in the case where two zero overlap values are present. Test case should yield the same results as the UNonOrthogonalStateBasis.
 */
BOOST_AUTO_TEST_CASE(NOCI_from_unrestricted_two_zeros) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    // It was for H3, at 1.6Å internuclear distance for the STO-3G basis set.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.8897259886);  // H3, 1.6Å apart.

    // The unrestricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    const auto S = spinor_basis.overlap();

    // Initialize two non-orthogonal "generalized states". These are constructed from the unrestricted case.
    // Columns 2 and 3 are swapped to make sure the correct occupied orbitals are selected by the algorithm.
    GQCP::SquareMatrix<double> state_1 {6};
    // clang-format off
    state_1 << 0.301612395,  0.0        ,  0.0        ,  1.18334666 ,  0.0        ,  0.0        ,
               0.460055205, -0.996503166,  0.0        , -0.535359699,  0.0        ,  0.0        , 
               0.460055205,  0.996503166,  0.0        , -0.535359699,  0.0        ,  0.0        ,
               0.0        ,  0.0        ,  0.638858962,  0.0        , -1.04073944 ,  0.0        ,
               0.0        ,  0.0        ,  0.280666711,  0.0        ,  0.647678158, -0.996503166,
               0.0        ,  0.0        ,  0.280666711,  0.0        ,  0.647678158,  0.996503166;
    // clang-format on
    GQCP::SquareMatrix<double> state_2 {6};
    // clang-format off
    state_2 <<  0.30161239,  0.63885896,  0.0        ,  1.18334666,  0.0        ,  0.0        ,
                0.46005521,  0.28066671,  0.0        , -0.5353597 ,  0.0        ,  0.0        ,
                0.46005521,  0.28066671,  0.0        , -0.5353597 ,  0.0        ,  0.0        ,
                0.0       ,  0.0       ,  0.0        ,  0.0       , -1.04073944 ,  0.0        ,
                0.0       ,  0.0       , -0.996503166,  0.0       ,  0.647678158, -0.996503166,
                0.0       ,  0.0       ,  0.996503166,  0.0       ,  0.647678158,  0.996503166;
    // clang-format on

    // Transform the matrices to the correct transformation type.
    const auto basis_state_1 = GQCP::GTransformation<double> {state_1};
    const auto basis_state_2 = GQCP::GTransformation<double> {state_2};

    // Create a vector out of these two basis states.
    std::vector<GQCP::GTransformation<double>> basis_vector {basis_state_1, basis_state_2};

    // Create a non-orthogonal state basis, using the basis state vector, the overlap operator in AO basis and the number of occupied orbitals.
    const auto NOS_basis = GQCP::GNonOrthogonalStateBasis<double> {basis_vector, S, 3};

    // To evaluate the Hamiltonian in the non-orthogonal state basis, we need the second quantized Hamiltonian in AO basis.
    const auto sq_hamiltonian = spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));

    // Initialize the reference overlap matrix.
    GQCP::SquareMatrix<double> reference_overlap {2};
    // clang-format off
    reference_overlap <<  1.0,  0.0      ,
                          0.0,  0.0878833;
    // clang-format on

    // Initialize the reference core hamiltonian.
    GQCP::SquareMatrix<double> reference_core_hamiltonian {2};
    // clang-format off
    reference_core_hamiltonian << -4.44852,   0.0     ,
                                   0.0    ,  -0.337126;
    // clang-format on

    // Initialize the reference two electron contribution.
    GQCP::SquareMatrix<double> reference_two_electron {2};
    // clang-format off
    reference_two_electron <<  1.52501    ,  -0.0196749,
                              -0.0196749  ,   0.131393 ;
    // clang-format on
    BOOST_CHECK(NOS_basis.evaluateOverlapOperator().isApprox(reference_overlap, 1e-6));
    BOOST_CHECK(NOS_basis.evaluateOneElectronOperator(sq_hamiltonian.core()).isApprox(reference_core_hamiltonian, 1e-5));
    BOOST_CHECK(NOS_basis.evaluateTwoElectronOperator(sq_hamiltonian.twoElectron()).isApprox(reference_two_electron, 1e-5));
}