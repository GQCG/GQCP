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

#define BOOST_TEST_MODULE "ScalarBasis"

#include <boost/test/unit_test.hpp>

#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Molecule/Molecule.hpp"


/**
 *  Check if the number of basis functions for H2O//STO-3G is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(numberOfBasisFunctions) {

    auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::ScalarBasis<GQCP::GTOShell> basis {molecule, "STO-3G"};

    BOOST_CHECK_EQUAL(basis.numberOfBasisFunctions(), 7);
}


/**
 *  Check if we can succesfully initialize a ScalarBasis for NO+ as large interatomic distances.
 */
BOOST_AUTO_TEST_CASE(dissociated_scalar_basis) {

    // Test if we can succesfully initialize NO+ at long intra molecular distance
    const GQCP::Nucleus N {7, 3.5, 0, 0};
    const GQCP::Nucleus O {8, -3.5, 0, 0};
    const std::vector<GQCP::Nucleus> nuclei {N, O};
    const auto molecule = GQCP::Molecule(nuclei, +1);

    BOOST_CHECK_NO_THROW(GQCP::ScalarBasis<GQCP::GTOShell>(molecule, "STO-3G"));
}


/**
 *  Check if the basis functions for H2O//STO-3G are correctly implemented.
 * 
 *  The example we're taking is a scalar basis derived from H2O//STO-3G.
 * 
 *  The relevant STO-3G basisset info is as follows:
 * 
 *      EXPONENTS              COEFFICIENTS (S)       COEFFICIENTS (P)
 * H
 * S
 *      3.42525091             0.15432897
 *      0.62391373             0.53532814
 *      0.16885540             0.44463454
 *
 * O
 * S
 *    130.7093200              0.15432897
 *     23.8088610              0.53532814
 *      6.4436083              0.44463454
 * SP
 *      5.0331513             -0.09996723             0.15591627
 *      1.1695961              0.39951283             0.60768372
 *      0.3803890              0.70011547             0.39195739
 */
BOOST_AUTO_TEST_CASE(basisFunctions) {

    // Set up the restricted spin-orbital basis from the basisset specification.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};

    const auto basis_functions = scalar_basis.basisFunctions();


    // In this set of basis functions, the 1py-function on O has the fourth index.
    const auto bf = basis_functions[4];

    // Check the coefficients and Gaussian exponents according to the basis set specification.
    // Since the contraction coefficients correspond to normalized primitives, we'll have to convert them to match unnormalized primitives since CartesianGTO is an unnormalized Cartesian Gaussian.
    std::vector<double> ref_coefficients {0.15591627, 0.60768372, 0.39195739};
    const std::vector<double> ref_gaussian_exponents {5.0331513, 1.1695961, 0.3803890};
    for (size_t d = 0; d < 3; d++) {
        ref_coefficients[d] *= GQCP::CartesianGTO::calculateNormalizationFactor(ref_gaussian_exponents[d], GQCP::CartesianExponents(0, 1, 0));
    }
    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(std::abs(bf.coefficient(i) - ref_coefficients[i]) < 1.0e-07);
        BOOST_CHECK(std::abs(bf.function(i).gaussianExponent() - ref_gaussian_exponents[i]) < 1.0e-07);
    }
}


/**
 *  Check if the basis function selector generates the correct indices for the AOs.
 */
BOOST_AUTO_TEST_CASE(selector_CO_basis_functions) {

    // Create a test molecule and scalar basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const GQCP::Vector<double, 3> C_position {0.0, 0.0, -1.23701};  // In Bohr.
    const GQCP::Vector<double, 3> O_position {0.0, 0.0, 0.927667};  // In Bohr.

    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};


    // Create the reference indices for CO in an STO-3G basis.
    const std::vector<size_t> ref_C_indices {0, 1, 2, 3, 4};  // The indices of the AOs centered on 'C'.
    const std::vector<size_t> ref_O_indices {5, 6, 7, 8, 9};  // The indices of the AOs centered on 'O'.


    // Check if the basis function selector generates the basis function correct indices.
    using BasisFunction = GQCP::ScalarBasis<GQCP::GTOShell>::BasisFunction;
    const auto C_indices = scalar_basis.basisFunctionIndices(
        [&C_position](const BasisFunction& bf) {
            return bf.functions()[0].center().isApprox(C_position, 1.0e-04);
        });
    BOOST_TEST(C_indices == ref_C_indices, boost::test_tools::per_element());


    const auto O_indices = scalar_basis.basisFunctionIndices(
        [&O_position](const BasisFunction& bf) {
            return bf.functions()[0].center().isApprox(O_position, 1.0e-04);
        });
    BOOST_TEST(O_indices == ref_O_indices, boost::test_tools::per_element());
}


/**
 *  Check if the shell selector generates the correct indices for the AOs.
 */
BOOST_AUTO_TEST_CASE(selector_CO_basis_shells) {

    // Create a test molecule and scalar basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/CO_mulliken.xyz");
    const GQCP::ScalarBasis<GQCP::GTOShell> scalar_basis {molecule, "STO-3G"};


    // Create the reference indices for CO in an STO-3G basis.
    const std::vector<size_t> ref_C_indices {0, 1, 2, 3, 4};  // The indices of the AOs centered on 'C'.
    const std::vector<size_t> ref_O_indices {5, 6, 7, 8, 9};  // The indices of the AOs centered on 'O'.

    // Check if the shell selector generates the correct basis function indices.
    using Shell = GQCP::ScalarBasis<GQCP::GTOShell>::Shell;
    const auto C_indices = scalar_basis.basisFunctionIndices(
        [](const Shell& shell) {
            return shell.nucleus().element() == "C";
        });
    BOOST_TEST(C_indices == ref_C_indices, boost::test_tools::per_element());


    const auto O_indices = scalar_basis.basisFunctionIndices(
        [](const Shell& shell) {
            return shell.nucleus().element() == "O";
        });
    BOOST_TEST(O_indices == ref_O_indices, boost::test_tools::per_element());
}
