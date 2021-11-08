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

#define BOOST_TEST_MODULE "GLowdinPairingBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/BiorthogonalBasis/GLowdinPairingBasis.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test whether or not the constructor of a biorthonormal `GLowdinPairingBasis` works.
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
    GQCP::SquareMatrix<double> C_ket {8};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
              0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
              0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
              0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
              0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
             -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
             -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::GTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states. Check whether the procedure doesn't throw any exceptions.
    const auto lowdin_pairing_basis = GQCP::GLowdinPairingBasis<double>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons());
    BOOST_CHECK_NO_THROW(GQCP::GLowdinPairingBasis<double>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons()));

    // Check whether the saved values are correct. The reference data is taken from the calculations of @lelemmen and @johdvos.
    // Initialize the reference biorthogonal overlaps.
    GQCP::VectorX<double> overlap_reference {2};
    // clang-format off
    overlap_reference <<  0.92270529,
                          0.84911655;
    // clang-format on

    // Initialize a reference for the biorthogonal occupied bra expansion coefficients.
    GQCP::MatrixX<double> biorthogonal_bra_reference {8, 2};
    // clang-format off
    biorthogonal_bra_reference <<  0.14038635, -0.01786507,
                                   0.16993508, -0.00914819,
                                   0.34156389, -0.11455247,
                                   0.53468019, -0.18444736,
                                  -0.11455247, -0.34156389,
                                  -0.18444736, -0.53468019,
                                  -0.01786507, -0.14038635,
                                  -0.00914819, -0.16993508;
    // clang-format off

    // Initialize a reference for the biorthogonal occupied ket expansion coefficients.
    GQCP::MatrixX<double> biorthogonal_ket_reference {8, 2};
    // clang-format off
    biorthogonal_ket_reference <<  0.18680645,  -0.23037179,
                                   0.23464362,  -0.40130146,
                                   0.29390833,   0.03997422,
                                   0.44583031,   0.13178613,
                                   0.00576840,  -0.25907694,
                                   0.02506212,  -0.38646563,
                                  -0.18158893,  -0.18484417,
                                  -0.28286936,  -0.26445598;
    // clang-format off

    // Check the parameters instantiated by the constructor.
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalBraExpansion().isApprox(biorthogonal_bra_reference, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalKetExpansion().isApprox(biorthogonal_ket_reference, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalOverlaps().isApprox(overlap_reference, 1e-6));
}


/**
 *  Test whether or not the constructor of a complex biorthonormal `GLowdinPairingBasis` works.
 */
BOOST_AUTO_TEST_CASE(constructor_complex) {
    // This is a self implemented test case, with rotated H3 states in STO-3G.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.88973, 0);  // H3, 1 Angstrom apart.

    // The general spinor basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::GSpinorBasis<GQCP::complex, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap();

    // Initialize two non-orthogonal "Generalised states".
    // Initialize some non-orthogonal "generalised states".
    GQCP::SquareMatrix<GQCP::complex> C_bra {6};
    // clang-format off
    C_bra << 0.460055771  - 5.63405828e-17j,  1.35013939  + 0.0j           ,  0.996501939    - 1.22036291e-16j,  0.535358621 - 6.55625222e-17j,  1.35013939 + 0.0j            ,  0.0 + 0.0j,
             0.460055767  - 5.63405822e-17j,  1.35013943  + 0.0j           , -0.996501918    + 1.22036288e-16j,  0.535358665 - 6.55625275e-17j,  1.35013943 + 0.0j            ,  0.0 + 0.0j,
             0.301611955  - 3.69368116e-17j,  3.07322180  + 0.0j           , -2.63772351e-08 + 3.23027965e-24j, -1.18334547  + 1.44918025e-16j,  3.07322180 + 0.0j            ,  0.0 + 0.0j,
            -0.0956361988 + 0.0j           , -0.280666403 - 3.43717213e-17j, -0.207152401    + 0.0j           , -0.111290123 + 0.0j           , -0.280666403 - 3.43717213e-17j,  0.0 + 0.0j,
            -0.0956361979 + 0.0j           , -0.280666413 - 3.43717224e-17j,  0.207152397    + 0.0j           , -0.111290132 + 0.0j           , -0.280666413 - 3.43717224e-17j,  0.0 + 0.0j,
            -0.0626989655 + 0.0j           , -0.638860046 - 7.82377910e-17j,  5.48328845e-09 + 0.0j           ,  0.245993356 + 0.0j           , -0.638860046 - 7.82377910e-17j,  0.0 + 0.0j;
    // clang-format on
    GQCP::SquareMatrix<GQCP::complex> C_ket {6};
    // clang-format off
    C_ket <<  0.230027886 - 0.398419985j,  0.799802113 + 0.0j        ,  0.498250970    - 0.862995994j   ,  0.267679311 - 0.463634166j,  0.799802113 + 0.0j        , 0.0 + 0.0j,
              0.230027883 - 0.398419981j,  0.799802140 + 0.0j        , -0.498250959    + 0.862995976j   ,  0.267679332 - 0.463634204j,  0.799802140 + 0.0j        , 0.0 + 0.0j,
              0.150805978 - 0.261203615j,  1.82053003  + 0.0j        , -1.31886175e-08 + 2.28433557e-08j, -0.591672737 + 1.02480724j ,  1.82053003  + 0.0j        , 0.0 + 0.0j,
             -0.161442683 + 0.0j        , -0.140333202 - 0.243064235j, -0.349692268    + 0.0j           , -0.187867944 + 0.0j        , -0.140333202 - 0.243064235j, 0.0 + 0.0j,
             -0.161442681 + 0.0j        , -0.140333206 - 0.243064243j,  0.349692261    + 0.0j           , -0.187867959 + 0.0j        , -0.140333206 - 0.243064243j, 0.0 + 0.0j,
             -0.105841609 + 0.0j        , -0.319430023 - 0.553269029j,  9.25629424e-09 + 0.0j           ,  0.415259365 + 0.0j        , -0.319430023 - 0.553269029j, 0.0 + 0.0j;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<GQCP::complex> {C_bra};
    const auto ket_expansion = GQCP::GTransformation<GQCP::complex> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states. Check whether the procedure doesn't throw any exceptions.
    const auto lowdin_pairing_basis = GQCP::GLowdinPairingBasis<GQCP::complex>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons());
    BOOST_CHECK_NO_THROW(GQCP::GLowdinPairingBasis<GQCP::complex>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons()));

    // Check whether the saved values are correct. The reference data is taken from the calculations of @lelemmen and @johdvos.
    // Initialize the reference biorthogonal overlaps.
    GQCP::VectorX<double> overlap_reference {3};
    // clang-format off
    overlap_reference <<  15.31141866,
                           1.03839811,
                           0.08484072;
    // clang-format on

    // Initialize a reference for the biorthogonal occupied bra expansion coefficients.
    GQCP::MatrixX<GQCP::complex> biorthogonal_bra_reference {6, 3};
    // clang-format off
    biorthogonal_bra_reference <<  -0.780318723 + 1.17946932j ,  0.566535554    - 0.819788749j   , -0.102459322  + 0.154869571j ,
                                   -0.780318747 + 1.17946935j , -0.566535535    + 0.819788722j   , -0.102459329  + 0.154869583j ,
                                   -1.69531930  + 2.56251329j ,  0.0            + 0.0j           ,  0.170136306  - 0.257164857j ,
                                    0.162212325 - 0.245187580j, -0.117771171    + 0.170417338j   ,  0.0212992006 - 0.0321942212j,
                                    0.162212331 - 0.245187589j,  0.117771167    - 0.170417333j   ,  0.0212992019 - 0.0321942233j,
                                    0.352422258 - 0.532694177j,  0.0            - 0.0j           , -0.0353678632 + 0.0534593215j;
    // clang-format off

    // Initialize a reference for the biorthogonal occupied ket expansion coefficients.
    GQCP::MatrixX<GQCP::complex> biorthogonal_ket_reference {6, 3};
    // clang-format off
    biorthogonal_ket_reference <<  -0.452370436 + 0.783528578j,  0.515610447 - 0.852737940j  , -0.0905244125 + 0.156792882j,
                                   -0.452370452 + 0.783528606j, -0.515610429 + 0.852737911j  , -0.0905244192 + 0.156792894j,
                                   -0.910842054 + 1.57762471j ,  0.0         + 0.0j          ,  0.147280666  - 0.255097596j,
                                    0.317491492 - 0.0j        , -0.349620681 - 0.00707551711j,  0.0635336192 - 0.0j        ,
                                    0.317491503 + 0.0j        ,  0.349620670 + 0.00707551697j,  0.0635336239 - 0.0j        ,
                                    0.639265034 + 0.0j        ,  0.0         + 0.0j          , -0.103367405  - 0.0j        ;
    // clang-format off

    // Check the parameters instantiated by the constructor.
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalBraExpansion().isApprox(biorthogonal_bra_reference, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalKetExpansion().isApprox(biorthogonal_ket_reference, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalOverlaps().isApprox(overlap_reference, 1e-6));
}


/**
 *  Test whether the overlap metrics are working correctly in the `GLowdinPairingBasis`.
 */
BOOST_AUTO_TEST_CASE(overlap_metrics) {
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
    GQCP::SquareMatrix<double> C_ket {8};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
              0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
              0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
              0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
              0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
             -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
             -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::GTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states. Check whether the procedure doesn't throw any exceptions.
    const auto lowdin_pairing_basis = GQCP::GLowdinPairingBasis<double>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons());

    // Initialize an empty reference vector, since there are no zero overlap values in this case.
    std::vector<int> zero_indices_ref {};

    // Check the overlap metrics calculated by the `LowdinPairingBasis'.
    BOOST_CHECK(std::abs(0.7834843312747031 - lowdin_pairing_basis.totalOverlap()) < 1e-6);
    BOOST_CHECK(std::abs(0.7834843312747031 - lowdin_pairing_basis.reducedOverlap()) < 1e-6);
    BOOST_CHECK_EQUAL(lowdin_pairing_basis.numberOfZeroOverlaps(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(zero_indices_ref.begin(), zero_indices_ref.end(), lowdin_pairing_basis.zeroOverlapIndices().begin(), lowdin_pairing_basis.zeroOverlapIndices().end());
}


/**
 *  Test whether the density matrices are working correctly in the `GLowdinPairingBasis`.
 */
BOOST_AUTO_TEST_CASE(density_matrices) {
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
    GQCP::SquareMatrix<double> C_ket {8};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,  0.81289875, -0.77260907,  0.50167389, -0.44422385,
              0.36593356, -0.28669343, -0.84796858, -0.13503625, -0.62437698,  0.96771154, -0.55231929,  0.30317456,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,  0.81450567,  0.7709918 , -0.501289  , -0.44451308,
              0.36597032,  0.28670189,  0.847938  , -0.13501526, -0.62639487, -0.96647128,  0.5520554 ,  0.30349133,
              0.10076798, -0.23874662,  0.04823423,  0.17685836,  0.42013282, -0.48352714, -0.79642816,  0.8239984 ,
              0.16561668, -0.35007843,  0.19502141,  0.90182842, -0.55545195,  0.39170258,  0.56753639, -0.94408827,
             -0.10075937, -0.23872464,  0.0482368 , -0.17686313, -0.42104909, -0.4826058 , -0.79588057, -0.82460595,
             -0.16560552, -0.35003836,  0.19503259, -0.9018579 ,  0.55619574,  0.39048771,  0.56690551,  0.94451894;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::GTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::GTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states. Check whether the procedure doesn't throw any exceptions.
    const auto lowdin_pairing_basis = GQCP::GLowdinPairingBasis<double>(bra_expansion, ket_expansion, S, molecule.numberOfElectrons());

    // Initialize a reference co-density matrix. Data taken from the implementation of @lelemmen and @johdvos.
    GQCP::SquareMatrix<double> ref_co_density {8};
    // clang-format off
    ref_co_density <<  2.62250750e-02,  3.17449693e-02,  6.38063371e-02,  9.98817070e-02, -2.13991402e-02, -3.44559565e-02, -3.33731112e-03, -1.70894095e-03,
                       3.29407612e-02,  3.98741835e-02,  8.01457883e-02,  1.25459296e-01, -2.68790067e-02, -4.32793970e-02, -4.19192580e-03, -2.14656452e-03,
                       4.12607161e-02,  4.99453355e-02,  1.00388470e-01,  1.57146958e-01, -3.36679246e-02, -5.42106148e-02, -5.25069411e-03, -2.68872929e-03,
                       6.25884885e-02,  7.57622103e-02,  1.52279534e-01,  2.38376633e-01, -5.10709632e-02, -8.22322238e-02, -7.96479167e-03, -4.07854051e-03,
                       8.09805133e-04,  9.80254169e-04,  1.97027842e-03,  3.08425120e-03, -6.60784900e-04, -1.06396685e-03, -1.03052963e-04, -5.27704553e-05,
                       3.51837883e-03,  4.25893264e-03,  8.56031357e-03,  1.34002165e-02, -2.87092722e-03, -4.62264103e-03, -4.47736557e-04, -2.29273001e-04,
                      -2.54926060e-02, -3.08583291e-02, -6.20242195e-02, -9.70920007e-02,  2.08014600e-02,  3.34935981e-02,  3.24409968e-03,  1.66121005e-03,
                      -3.97109961e-02, -4.80694279e-02, -9.66179583e-02, -1.51244642e-01,  3.24033839e-02,  5.21745068e-02,  5.05348217e-03,  2.58774273e-03;
    // clang-format on

    // Initialize a zero matrix as reference for the zero overlap co-density matrix.
    const GQCP::MatrixX<double> zero_overlap_co_dens_ref = GQCP::MatrixX<double>::Zero(8, 8);

    // Initialize a reference weighted co-density matrix. Data taken from the implementation of @lelemmen and @johdvos.
    GQCP::SquareMatrix<double> ref_weighted_co_density {8};
    // clang-format off
    ref_weighted_co_density << 0.03326887,  0.03688621,  0.10023033,  0.15829074,  0.06947715,  0.10772049,  0.03447101,  0.04425258,
                               0.04414342,  0.04753796,  0.14099828,  0.22314075,  0.13229605,  0.20579061,  0.06180499,  0.07798674,
                               0.04387607,  0.05369857,  0.10340514,  0.1616278 , -0.05256822, -0.08392318, -0.01229957, -0.01091407,
                               0.06505878,  0.08068895,  0.14725696,  0.22971843, -0.10836118, -0.17210519, -0.0304205 , -0.03079477,
                               0.00632852,  0.00385361,  0.03708684,  0.05962   ,  0.10349962,  0.16198506,  0.04272209,  0.0517923 ,
                               0.0119442 ,  0.0087794 ,  0.06141465,  0.09847184,  0.15234745,  0.23834366,  0.06340999,  0.07709552,
                              -0.02373907, -0.03145185, -0.04228304, -0.06507302,  0.09689903,  0.15269387,  0.03407656,  0.03879353,
                              -0.03747352, -0.049247  , -0.06903443, -0.10646852,  0.14149733,  0.22307043,  0.04919991,  0.05573053;
    // clang-format on

    // Initialize a reference weighted co-density matrix. Data taken from the implementation of @lelemmen and @johdvos.
    GQCP::SquareMatrix<double> ref_transition_1DM {8};
    // clang-format off
    ref_transition_1DM << 0.02606564,  0.03458568,  0.03437621,  0.05097253,  0.00495829,  0.00935809, -0.01859919, -0.02935992,
                          0.02889977,  0.03724525,  0.04207199,  0.06321853,  0.00301924,  0.00687852, -0.02464203, -0.03858425,
                          0.07852889,  0.11046995,  0.08101631,  0.11537352,  0.02905695,  0.04811741, -0.0331281 , -0.05408739,
                          0.12401831,  0.17482728,  0.12663285,  0.17998079,  0.04671134,  0.07715115, -0.05098369, -0.08341642,
                          0.05443426,  0.10365189, -0.04118638, -0.08489929,  0.08109033,  0.11936184,  0.07591887,  0.11086094,
                          0.08439732,  0.16123372, -0.0657525 , -0.13484172,  0.12691275,  0.18673852,  0.11963325,  0.17477219,
                          0.0270075 ,  0.04842324, -0.00963652, -0.02383398,  0.03347209,  0.04968073,  0.02669845,  0.03854736,
                          0.03467121,  0.06110139, -0.008551  , -0.02412722,  0.04057846,  0.06040313,  0.03039413,  0.04366399;
    // clang-format on

    // Compare the calculated and reference co-density matrices.
    BOOST_CHECK(lowdin_pairing_basis.coDensity(0).matrix().isApprox(ref_co_density, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.zeroOverlapCoDensity().matrix().isApprox(zero_overlap_co_dens_ref, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.weightedCoDensity().matrix().isApprox(ref_weighted_co_density, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.transition1DM().matrix().isApprox(ref_transition_1DM, 1e-6));

    // The co-density matrix sum yields the same result as the regular weighted co-density matrix as there are no zero overlaps in this calculation.
    BOOST_CHECK(lowdin_pairing_basis.coDensitySum().matrix().isApprox(ref_weighted_co_density, 1e-6));
}
