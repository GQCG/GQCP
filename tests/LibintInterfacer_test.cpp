// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "LibintInterfacer"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "LibintInterfacer.hpp"

#include "utilities/linalg.hpp"


/*
 *  UTILITIES
 */

/**
 *  A functor to compare two doubles with respect to a tolerance
 */
struct approx {
public:
    double tolerance;


public:
    approx(double tolerance = 1.0e-12) : tolerance(tolerance) {}

    bool operator()(double lhs, double rhs) const {
        return std::abs(lhs - rhs) < tolerance;
    }
};




/*
 *  UNIT TESTS
 */

BOOST_AUTO_TEST_CASE ( atoms_to_libint ) {

    std::vector<GQCP::Atom> GQCP_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    std::vector<libint2::Atom> ref_libint_atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };


    // Use the Libint interface to obtain a std::vector<libint2::Atom> from the GQCP ones
    auto test_libint_atoms = GQCP::LibintInterfacer::get().interface(GQCP_atoms);


    /**
     *  Implement a function object to compare (libint_atom) == (libint_atom)
     */
    struct LibintAtomEqual {
        double tolerance;
        explicit LibintAtomEqual(double tolerance) : tolerance (tolerance) {};

        bool operator()(const libint2::Atom& lhs, const libint2::Atom& rhs) {
            return (lhs.atomic_number == rhs.atomic_number) &&
                   (std::abs(lhs.x - rhs.x) < tolerance) &&
                   (std::abs(lhs.y - rhs.y) < tolerance) &&
                   (std::abs(lhs.z - rhs.z) < tolerance);
        }
    };


    // Check if the interfacing between the Atom types works
    BOOST_CHECK((ref_libint_atoms.size() == test_libint_atoms.size()) &&
                std::equal(ref_libint_atoms.begin(), ref_libint_atoms.end(), test_libint_atoms.begin(), LibintAtomEqual(1.0e-08)));
}


BOOST_AUTO_TEST_CASE ( Shell_to_libint ) {

    // Create a GQCP::Shell
    GQCP::Atom h (1,  0.0, 0.0, 0.0);
    GQCP::Shell shell (0, h, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, false);  // a Cartesian s-type shell on the H-atom

    // Create the reference libint2::Shell (indirectly, because we want to refrain from renorm() being called
    libint2::Shell ref_libint_shell {};
    ref_libint_shell.alpha = {3.42525091, 0.62391373, 0.16885540};
    ref_libint_shell.contr = {{0, false, {0.15432897, 0.53532814, 0.44463454}}};
    ref_libint_shell.O = {0.0, 0.0, 0.0};


    // Check if the interfacing from a GQCP::Shell to a libint2::Shell works
    // Note that this function also tests LibintInterfacer::undo_renorm()
    auto libint_shell = GQCP::LibintInterfacer::get().interface(shell);

    BOOST_CHECK(std::equal(ref_libint_shell.alpha.begin(), ref_libint_shell.alpha.end(), libint_shell.alpha.begin(), approx()));  // check Gaussian exponents
    BOOST_CHECK(std::equal(ref_libint_shell.contr[0].coeff.begin(), ref_libint_shell.contr[0].coeff.end(), libint_shell.contr[0].coeff.begin(), approx()));  // check contraction coefficients
    BOOST_CHECK(ref_libint_shell.contr[0].pure == libint_shell.contr[0].pure);  // check if both shells have the same type
    BOOST_CHECK(std::equal(ref_libint_shell.O.begin(), ref_libint_shell.O.end(), libint_shell.O.begin(), approx()));  // check the center
}


BOOST_AUTO_TEST_CASE ( ShellSet_to_BasisSet ) {

    // Since we can't really create an initialized libint2::BasisSet except when using a molecule and a basisset name, we can only test if the shells are correctly placed and if the initialization 'hack' does what is expected


    // Make a test ShellSet
    GQCP::Shell s (0, GQCP::Atom(), {1.0, 2.0}, {0.5, -0.5}, true);
    GQCP::Shell d (2, GQCP::Atom(), {4.0, 8.0}, {1.5, -1.5}, false);
    GQCP::ShellSet shellset {s, d};


    // Make the reference BasisSet
    libint2::BasisSet ref_libint_basisset {};
    ref_libint_basisset.push_back(libint2::Shell());
    ref_libint_basisset[0].alpha = {1.0, 2.0};
    ref_libint_basisset[0].contr = {{0, true, {0.5, -0.5}}};
    ref_libint_basisset[0].O = {0.0, 0.0, 0.0};

    ref_libint_basisset.push_back(libint2::Shell());
    ref_libint_basisset[1].alpha = {4.0, 8.0};
    ref_libint_basisset[1].contr = {{2, false, {1.5, -1.5}}};
    ref_libint_basisset[1].O = {0.0, 0.0, 0.0};


    // Check if the interfacing works
    auto libint_basisset = GQCP::LibintInterfacer::get().interface(shellset);

    for (size_t i = 0; i < ref_libint_basisset.size(); i++) {
        BOOST_CHECK(std::equal(ref_libint_basisset[i].alpha.begin(), ref_libint_basisset[i].alpha.end(), libint_basisset[i].alpha.begin(), approx()));  // check Gaussian exponents
        BOOST_CHECK(std::equal(ref_libint_basisset[i].contr[0].coeff.begin(), ref_libint_basisset[i].contr[0].coeff.end(), libint_basisset[i].contr[0].coeff.begin(), approx()));  // check contraction coefficients
        BOOST_CHECK(ref_libint_basisset[i].contr[0].pure == libint_basisset[i].contr[0].pure);  // check if both shells have the same type
        BOOST_CHECK(std::equal(ref_libint_basisset[i].O.begin(), ref_libint_basisset[i].O.end(), libint_basisset[i].O.begin(), approx()));  // check the center
    }

    // Check if the initialization 'hack' did something
    BOOST_CHECK(libint_basisset.nbf() != -1);
    BOOST_CHECK(libint_basisset.max_nprim() != 0);
    BOOST_CHECK(libint_basisset.max_l() != -1);
}


BOOST_AUTO_TEST_CASE ( BasisSet_to_ShellSet ) {

    // Note that this function also tests std::vector<Shell> LibintInterfacer::interface(const libint2::Shell& libint_shell, const std::vector<Atom>& atoms) const;

    
    // Create a reference STO-3G shell set on (a weird geometry of) H2O
    GQCP::Atom h1 (1,  0.0, 0.0, 0.0);
    GQCP::Atom o  (8,  0.0, 0.0, 1.0);
    GQCP::Atom h2 (1,  0.0, 0.0, 2.0);
    std::vector<GQCP::Atom> atoms {h1, o, h2};

    GQCP::ShellSet ref_shellset {
        GQCP::Shell(0, h1, {  3.42525091,  0.62391373, 0.16885540}, { 0.15432897, 0.53532814, 0.44463454}),
        GQCP::Shell(0, o,  {130.7093200,  23.8088610,  6.4436083},  { 0.15432897, 0.53532814, 0.44463454}),
        GQCP::Shell(0, o,  {  5.0331513,   1.1695961,  0.3803890},  {-0.09996723, 0.39951283, 0.70011547}),
        GQCP::Shell(1, o,  {  5.0331513,   1.1695961,  0.3803890},  { 0.15591627, 0.60768372, 0.39195739}),
        GQCP::Shell(0, h2, {  3.42525091,  0.62391373, 0.16885540}, { 0.15432897, 0.53532814, 0.44463454})
    };

    GQCP::Molecule h2o ({h1, o, h2});


    // Construct the corresponding libint2::BasisSet
    auto libint_atoms = GQCP::LibintInterfacer::get().interface(atoms);
    libint2::BasisSet libint_basisset ("STO-3G", libint_atoms);
    auto shellset = GQCP::LibintInterfacer::get().interface(libint_basisset, atoms);

    BOOST_CHECK(ref_shellset == shellset);
}
