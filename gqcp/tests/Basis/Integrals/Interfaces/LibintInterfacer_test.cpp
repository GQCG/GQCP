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

#define BOOST_TEST_MODULE "LibintInterfacer"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Mathematical/Representation/Matrix.hpp"


/*
 *  MARK: Utilities for the unit tests
 */

/**
 *  A functor to compare two doubles with respect to a tolerance.
 */
struct approx {
public:
    double tolerance;


public:
    approx(double tolerance = 1.0e-12) :
        tolerance(tolerance) {}

    bool operator()(double lhs, double rhs) const {
        return std::abs(lhs - rhs) < tolerance;
    }
};


/*
 *  MARK: The actual unit tests
 */

/**
 *  Check if the interfacing between libint2::Atom and GQCP::Nucleus works as expected.
 */
BOOST_AUTO_TEST_CASE(nuclei_to_libint) {

    // Initialize a set of GQCP nuclei and libint atoms.
    const std::vector<GQCP::Nucleus> GQCP_nuclei {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};

    const std::vector<libint2::Atom> ref_libint_atoms {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}};


    // Use the Libint interface to obtain a std::vector<libint2::Atom> from the GQCP nuclei.
    const auto test_libint_atoms = GQCP::LibintInterfacer::get().interface(GQCP_nuclei);


    /**
     *  A functor that compares compare (libint_atom) == (libint_atom)
     */
    struct LibintAtomEqual {
        double tolerance;
        explicit LibintAtomEqual(double tolerance) :
            tolerance(tolerance) {};

        bool operator()(const libint2::Atom& lhs, const libint2::Atom& rhs) {
            return (lhs.atomic_number == rhs.atomic_number) &&
                   (std::abs(lhs.x - rhs.x) < tolerance) &&
                   (std::abs(lhs.y - rhs.y) < tolerance) &&
                   (std::abs(lhs.z - rhs.z) < tolerance);
        }
    };


    // Check if the interfacing between the atom types behaves correctly.
    BOOST_CHECK((ref_libint_atoms.size() == test_libint_atoms.size()) &&
                std::equal(ref_libint_atoms.begin(), ref_libint_atoms.end(), test_libint_atoms.begin(), LibintAtomEqual(1.0e-08)));
}


/**
 *  Check if the interfacing between a GQCP::GTOShell and a libint2::Shell works as expected.
 */
BOOST_AUTO_TEST_CASE(Shell_to_libint) {

    // Create a GQCP::GTOShell. A Cartesian s-type shell on an H atom.
    const GQCP::Nucleus h {1, 0.0, 0.0, 0.0};
    const GQCP::GTOShell shell {0, h, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, false};

    // Create the reference libint2::GTOShell (indirectly, because we want to refrain from renorm() being called.
    libint2::Shell ref_libint_shell {};
    ref_libint_shell.alpha = {3.42525091, 0.62391373, 0.16885540};
    ref_libint_shell.contr = {{0, false, {0.15432897, 0.53532814, 0.44463454}}};
    ref_libint_shell.O = {0.0, 0.0, 0.0};


    // Check if the interfacing from a GQCP::GTOShell to a libint2::GTOShell works.
    // Note that this function also tests LibintInterfacer::undoRenorm().
    const auto libint_shell = GQCP::LibintInterfacer::get().interface(shell);

    BOOST_CHECK(std::equal(ref_libint_shell.alpha.begin(), ref_libint_shell.alpha.end(),
                           libint_shell.alpha.begin(), approx()));  // Check the Gaussian exponents.
    BOOST_CHECK(std::equal(ref_libint_shell.contr[0].coeff.begin(), ref_libint_shell.contr[0].coeff.end(),
                           libint_shell.contr[0].coeff.begin(), approx()));     // Check the contraction coefficients.
    BOOST_CHECK(ref_libint_shell.contr[0].pure == libint_shell.contr[0].pure);  // Check if both shells have the same type.
    BOOST_CHECK(std::equal(ref_libint_shell.O.begin(), ref_libint_shell.O.end(),
                           libint_shell.O.begin(), approx()));  // Check the center.
}


/**
 *  Check if the interfacing between a libint2::Shell and a GQCP::GTOShell works as expected.
 */
BOOST_AUTO_TEST_CASE(libint_Shell_to_Shell) {

    // Make some libint Shells.
    libint2::Shell libint_s_shell {};
    libint_s_shell.alpha = {1.0, 2.0};
    libint_s_shell.contr = {{0, true, {0.5, -0.5}}};
    libint_s_shell.O = {0.0, 0.0, 0.0};

    libint2::Shell libint_d_shell {};
    libint_d_shell.alpha = {4.0, 8.0};
    libint_d_shell.contr = {{2, false, {1.5, -1.5}}};
    libint_d_shell.O = {0.0, 0.0, 0.0};


    // Create the reference GQCP Shells.
    const GQCP::GTOShell ref_s_shell {0, GQCP::Nucleus(), {1.0, 2.0}, {0.5, -0.5}, true};
    const GQCP::GTOShell ref_d_shell {2, GQCP::Nucleus(), {4.0, 8.0}, {1.5, -1.5}, false};


    // Convert the libint2::Shells to GQCP::GTOShells and check if they match the references.
    const bool undo_renorm = false;
    const auto s_shell = GQCP::LibintInterfacer::get().interface(libint_s_shell, {GQCP::Nucleus()}, undo_renorm)[0];
    const auto d_shell = GQCP::LibintInterfacer::get().interface(libint_d_shell, {GQCP::Nucleus()}, undo_renorm)[0];

    BOOST_CHECK(ref_s_shell == s_shell);
    BOOST_CHECK(ref_d_shell == d_shell);
}


/**
 *  Check if the interfacing between a GQCP::ShellSet and a libint2::BasisSet works as expected.
 */
BOOST_AUTO_TEST_CASE(ShellSet_to_BasisSet) {

    // Initialize a GQCP::ShellSet.
    const GQCP::GTOShell s {0, GQCP::Nucleus(), {1.0, 2.0}, {0.5, -0.5}, true};   // An s-type shell.
    const GQCP::GTOShell d {2, GQCP::Nucleus(), {4.0, 8.0}, {1.5, -1.5}, false};  // A d-type shell.
    GQCP::ShellSet<GQCP::GTOShell> shellset {s, d};


    // Initialize the reference libint2::BasisSet.
    libint2::BasisSet ref_libint_basisset {};
    ref_libint_basisset.push_back(libint2::Shell());
    ref_libint_basisset[0].alpha = {1.0, 2.0};
    ref_libint_basisset[0].contr = {{0, true, {0.5, -0.5}}};
    ref_libint_basisset[0].O = {0.0, 0.0, 0.0};

    ref_libint_basisset.push_back(libint2::Shell());
    ref_libint_basisset[1].alpha = {4.0, 8.0};
    ref_libint_basisset[1].contr = {{2, false, {1.5, -1.5}}};
    ref_libint_basisset[1].O = {0.0, 0.0, 0.0};


    // Check if the interfacing works.
    auto libint_basisset = GQCP::LibintInterfacer::get().interface(shellset);

    for (size_t i = 0; i < ref_libint_basisset.size(); i++) {
        BOOST_CHECK(std::equal(ref_libint_basisset[i].alpha.begin(), ref_libint_basisset[i].alpha.end(),
                               libint_basisset[i].alpha.begin(), approx()));  // Check the Gaussian exponents.
        BOOST_CHECK(std::equal(ref_libint_basisset[i].contr[0].coeff.begin(), ref_libint_basisset[i].contr[0].coeff.end(),
                               libint_basisset[i].contr[0].coeff.begin(), approx()));           // check contraction coefficients
        BOOST_CHECK(ref_libint_basisset[i].contr[0].pure == libint_basisset[i].contr[0].pure);  // check if both shells have the same type
        BOOST_CHECK(std::equal(ref_libint_basisset[i].O.begin(), ref_libint_basisset[i].O.end(),
                               libint_basisset[i].O.begin(), approx()));  // check the center
    }
}


/**
 *  Check if a libint2::BasisSet is correctly converted to a vector of GQCP::GTOShells.
 */
BOOST_AUTO_TEST_CASE(BasisSet_to_ShellSet) {

    // Create a reference STO-3G shell set on (a weird geometry of) H2O.
    const GQCP::Nucleus h1 {1, 0.0, 0.0, 0.0};
    const GQCP::Nucleus o {8, 0.0, 0.0, 1.0};
    const GQCP::Nucleus h2 {1, 0.0, 0.0, 2.0};
    const std::vector<GQCP::Nucleus> nuclei {h1, o, h2};
    const bool pure = false;  // We'll use an STO-3G basis set, which represents Cartesian shells.

    const GQCP::ShellSet<GQCP::GTOShell> ref_shellset {
        GQCP::GTOShell(0, h1, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, pure),
        GQCP::GTOShell(0, o, {130.7093200, 23.8088610, 6.4436083}, {0.15432897, 0.53532814, 0.44463454}, pure),
        GQCP::GTOShell(0, o, {5.0331513, 1.1695961, 0.3803890}, {-0.09996723, 0.39951283, 0.70011547}, pure),
        GQCP::GTOShell(1, o, {5.0331513, 1.1695961, 0.3803890}, {0.15591627, 0.60768372, 0.39195739}, pure),
        GQCP::GTOShell(0, h2, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, pure)};


    // Create a libint2::BasisSet and check if it is correctly converted to a vector of GQCP::GTOShells.
    const auto libint_atoms = GQCP::LibintInterfacer::get().interface(nuclei);
    const libint2::BasisSet libint_basisset {"STO-3G", libint_atoms};

    const auto shells = GQCP::LibintInterfacer::get().interface(libint_basisset, nuclei);
    BOOST_CHECK(ref_shellset.asVector() == shells);
}
