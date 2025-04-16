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

#define BOOST_TEST_MODULE "GTOEvalAtPoint"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"

namespace tt = boost::test_tools;


/**
 *  Check if the GTO evaluation at a random point in space is identical to PySCF results.
 */
BOOST_AUTO_TEST_CASE(eval_gto) {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 1.0, 0.0, 0.0};
    const GQCP::Nucleus h2 {1, 2.0, 0.0, 0.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    // init basis
    const auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> {molecule, "STO-3G"};
    // gather AOs
    // each spatial orbital ("spatial MO") is expanded in the AO basis set functions.
    const auto AOs = spin_orbital_basis.spatialOrbitals()[0].functions();
    const int nao = AOs.size();

    // init random point in space
    const GQCP::Vector<double, 3> r = {1.3, -0.9, 3.7};

    // init vector
    std::vector<double> AO_vals;
    AO_vals.reserve(nao);
    // loop through AOs
    for (size_t i = 0; i < nao; i++) {
        // gather value at r, put it in vector
        AO_vals.push_back(AOs[i](r));
    }

    // compare with hardcoded reference data
    const std::vector<double> AO_vals_ref = {0.0054346, 0.0, 0.000939672, 0.000194685, -0.000584055, 0.002401113, 0.006665074};
    for (size_t i = 0; i < nao; i++) {
        BOOST_TEST(AO_vals[i] == AO_vals_ref[i], tt::tolerance(1e-6));
    }
}


/**
 *  Check if the London GTO evaluation at a random point in space is identical to PySCF results.
 */
 BOOST_AUTO_TEST_CASE(eval_london_gto) {

    // Create an STO-3G basisset on (a toy geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 1.0, 0.0, 0.0};
    const GQCP::Nucleus h2 {1, 2.0, 0.0, 0.0};
    const GQCP::Molecule molecule {{h1, o, h2}};

    // zero-field because PySCF does not have field-dependent AOs
    const auto B = GQCP::HomogeneousMagneticField {{0.0, 0.0, 0.0}};  // Gauge origin at cartesian origin.
    // init basis
    auto spin_orbital_basis = GQCP::RSpinOrbitalBasis<GQCP::complex, GQCP::LondonGTOShell> {molecule, "STO-3G", B};
    // gather AOs
    // each spatial orbital ("spatial MO") is expanded in the AO basis set functions.
    // const auto AOs = spin_orbital_basis.spatialOrbitals()[0].functions();
    const auto AOs = spin_orbital_basis.spatialOrbitals()[0].functions();
    const int nao = AOs.size();

    // init random point in space
    const GQCP::Vector<double, 3> r = {1.3, -0.9, 3.7};

    // init vector
    std::vector<std::complex<double>> AO_vals;
    AO_vals.reserve(nao);
    // loop through AOs
    for (size_t i = 0; i < nao; i++) {
        // gather value at r, put it in vector
        AO_vals.push_back(AOs[i](r));
    }

    // compare with hardcoded reference data
    const std::vector<std::complex<double>> AO_vals_ref = {
        {0.0054346, 0.0}, {0.0, 0.0}, {0.000939672, 0.0}, {0.000194685, 0.0}, {-0.000584055, 0.0}, {0.002401113, 0.0}, {0.006665074, 0.0}
    };
    for (size_t i = 0; i < nao; ++i) {
        BOOST_TEST(AO_vals[i].real() == AO_vals_ref[i].real(), tt::tolerance(1e-6));
        BOOST_TEST(AO_vals[i].imag() == AO_vals_ref[i].imag(), tt::tolerance(1e-6));
    }
}
