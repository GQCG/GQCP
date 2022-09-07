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

#define BOOST_TEST_MODULE "RLowdinPairingBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/BiorthogonalBasis/RLowdinPairingBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Test whether or not the constructor of a biorthonormal `RLowdinPairingBasis` works.
 */
BOOST_AUTO_TEST_CASE(constructor) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S_res = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "Restricted states".
    GQCP::SquareMatrix<double> C_bra {4};
    // clang-format off
    C_bra << -0.07443693,  0.12036042, -0.13557067,  0.15517005,
             -0.07874922,  0.15086478, -0.68085546,  0.77423311,
             -0.24580188,  0.26338108,  0.09556297, -0.12178159,
             -0.38944259,  0.4101685 ,  0.45214166, -0.58335985;
    // clang-format on
    GQCP::SquareMatrix<double> C_ket {4};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,
              0.36593356, -0.28669343, -0.84796858, -0.13503625,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,
              0.36597032,  0.28670189,  0.847938  , -0.13501526;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::RTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::RTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states. Check whether the procedure doesn't throw any exceptions.
    const auto lowdin_pairing_basis = GQCP::RLowdinPairingBasis<double>(bra_expansion, ket_expansion, S_res, molecule.numberOfElectronPairs());
    BOOST_CHECK_NO_THROW(GQCP::RLowdinPairingBasis<double>(bra_expansion, ket_expansion, S_res, molecule.numberOfElectronPairs()));

    // Initialize a reference for the biorthogonal occupied bra expansion coefficients.
    // Reference code results in two completely positive sets of expansion coefficients. However, as both phases change, this has no effect on the resulting calculations (as shown by  all subsequent tests).
    GQCP::MatrixX<double> biorthogonal_bra_reference {4, 1};
    // clang-format off
    biorthogonal_bra_reference <<  -0.07443693,
                                   -0.07874922,
                                   -0.24580188,
                                   -0.38944259;
    // clang-format off

    // Initialize a reference for the biorthogonal occupied ket expansion coefficients.
    GQCP::MatrixX<double> biorthogonal_ket_reference {4, 1};
    // clang-format off
    biorthogonal_ket_reference <<  0.25851329,
                                   0.36593356,
                                   0.25853403,
                                   0.36597032;
    // clang-format off

    // Check the parameters instantiated by the constructor.
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalBraExpansion().isApprox(biorthogonal_bra_reference, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.biorthogonalKetExpansion().isApprox(biorthogonal_ket_reference, 1e-6));

    // Due to the dimensions of the problem, there is only one biorthogonal overlap calculated. Since it is stored in a vector, we can compare it to the reference value as follows.
    BOOST_CHECK(std::abs(-0.591737 - lowdin_pairing_basis.biorthogonalOverlaps()[0]) < 1e-6);
}


/**
 *  Test whether the overlap metrics are working correctly in the `RLowdinPairingBasis`.
 */
BOOST_AUTO_TEST_CASE(overlap_metrics) {
    // Reference data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S_res = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "Restricted states".
    GQCP::SquareMatrix<double> C_bra {4};
    // clang-format off
    C_bra << -0.07443693,  0.12036042, -0.13557067,  0.15517005,
             -0.07874922,  0.15086478, -0.68085546,  0.77423311,
             -0.24580188,  0.26338108,  0.09556297, -0.12178159,
             -0.38944259,  0.4101685 ,  0.45214166, -0.58335985;
    // clang-format on
    GQCP::SquareMatrix<double> C_ket {4};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,
              0.36593356, -0.28669343, -0.84796858, -0.13503625,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,
              0.36597032,  0.28670189,  0.847938  , -0.13501526;
    // clang-format on
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::RTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::RTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states.
    const auto lowdin_pairing_basis = GQCP::RLowdinPairingBasis<double>(bra_expansion, ket_expansion, S_res, molecule.numberOfElectronPairs());

    // Initialize an empty reference vector, since there are no zero overlap values in this case.
    std::vector<int> zero_indices_ref {};

    // Check the overlap metrics calculated by the LowdinPairingBasis.
    BOOST_CHECK(std::abs(-0.591737 - lowdin_pairing_basis.totalOverlap()) < 1e-6);
    BOOST_CHECK(std::abs(-0.591737 - lowdin_pairing_basis.reducedOverlap()) < 1e-6);
    BOOST_CHECK_EQUAL(lowdin_pairing_basis.numberOfZeroOverlaps(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(zero_indices_ref.begin(), zero_indices_ref.end(), lowdin_pairing_basis.zeroOverlapIndices().begin(), lowdin_pairing_basis.zeroOverlapIndices().end());
}


/**
 *  Test whether the density matrices are working correctly in the `RLowdinPairingBasis`.
 */
BOOST_AUTO_TEST_CASE(density_matrices) {
    // It was for H2, at 2.5au internuclear distance for the 6-31G basis set.
    const auto molecule = GQCP::Molecule::HChain(2, 2.5, 0);  // H2, 2.5 bohr apart.

    // The restricted spin orbital basis is also needed, as we require the overlap operator in AO basis.
    const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "6-31G"};
    const auto S_res = spin_orbital_basis.overlap();

    // Initialize two non-orthogonal "restricted states".
    GQCP::SquareMatrix<double> C_bra {4};
    // clang-format off
    C_bra << -0.07443693,  0.12036042, -0.13557067,  0.15517005,
             -0.07874922,  0.15086478, -0.68085546,  0.77423311,
             -0.24580188,  0.26338108,  0.09556297, -0.12178159,
             -0.38944259,  0.4101685 ,  0.45214166, -0.58335985;
    // clang-format on
    GQCP::SquareMatrix<double> C_ket {4};
    // clang-format off
    C_ket <<  0.25851329, -0.14539151, -0.17177142, -0.01126487,
              0.36593356, -0.28669343, -0.84796858, -0.13503625,
              0.25853403,  0.14539669,  0.17176599, -0.01126146,
              0.36597032,  0.28670189,  0.847938  , -0.13501526;
    // clang-format on
    // Transform the matrices to the correct transformation type.
    const auto bra_expansion = GQCP::RTransformation<double> {C_bra};
    const auto ket_expansion = GQCP::RTransformation<double> {C_ket};

    // We can now initialize the biorthogonal Löwdin Pairing basis of these two states.
    const auto lowdin_pairing_basis = GQCP::RLowdinPairingBasis<double>(bra_expansion, ket_expansion, S_res, molecule.numberOfElectronPairs());

    // Initialize a reference co-density matrix. Data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    GQCP::SquareMatrix<double> ref_co_density {4};
    // clang-format off
    ref_co_density <<  -0.01924294, -0.02035772, -0.06354305, -0.10067609,
                       -0.02723897, -0.02881698, -0.08994716, -0.14251011,
                       -0.01924448, -0.02035935, -0.06354815, -0.10068416,
                       -0.02724171, -0.02881988, -0.08995619, -0.14252443;
    // clang-format on

    // Initialize a zero matrix as reference for the zero overlap co-density matrix.
    const GQCP::MatrixX<double> zero_overlap_co_dens_ref = GQCP::MatrixX<double>::Zero(4, 4);

    // Initialize a reference weighted co-density matrix. Data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    GQCP::SquareMatrix<double> ref_weighted_co_density {4};
    // clang-format off
    ref_weighted_co_density << 0.03251942, 0.03440334, 0.10738399, 0.17013661,
                               0.04603225, 0.04869899, 0.15200536, 0.24083364,
                               0.03252203, 0.0344061 , 0.10739261, 0.17015026,
                               0.04603687, 0.04870388, 0.15202063, 0.24085783;
    // clang-format on

    // Initialize a reference weighted co-density matrix. Data taken from the implementation of H. Burton (https://github.com/hgaburton/libgnme).
    GQCP::SquareMatrix<double> ref_transition_1DM {4};
    // clang-format off
    ref_transition_1DM << -0.01924294, -0.02723897, -0.01924448, -0.02724171,
                          -0.02035772, -0.02881698, -0.02035935, -0.02881988,
                          -0.06354305, -0.08994716, -0.06354815, -0.08995619,
                          -0.10067609, -0.14251011, -0.10068416, -0.14252443;
    // clang-format on

    // Compare the calculated and reference co-density matrices.
    BOOST_CHECK(lowdin_pairing_basis.coDensity(0).matrix().isApprox(ref_co_density, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.zeroOverlapCoDensity().matrix().isApprox(zero_overlap_co_dens_ref, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.weightedCoDensity().matrix().isApprox(ref_weighted_co_density, 1e-6));
    BOOST_CHECK(lowdin_pairing_basis.transition1DM().matrix().isApprox(ref_transition_1DM, 1e-6));

    // The co-density matrix sum yields the same result as the regular weighted co-density matrix as there are no zero overlaps in this calculation.
    BOOST_CHECK(lowdin_pairing_basis.coDensitySum().matrix().isApprox(ref_weighted_co_density, 1e-6));
}
