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


BOOST_AUTO_TEST_CASE ( shell_to_libint ) {

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


    // Create a GQCP::Shell and its reference
    GQCP::Atom h (1,  0.0, 0.0, 0.0);
    GQCP::Shell shell (0, h, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}, false);  // a Cartesian s-type shell on the H-atom

    libint2::Shell ref_libint_shell {
        {3.42525091, 0.62391373, 0.16885540},
        {
            {0, false, {0.15432897, 0.53532814, 0.44463454}}
        },
        {{0.0, 0.0, 0.0}}
    };


    // Check if the interfacing from a GQCP::Shell to a libint2::Shell works
    auto libint_shell = GQCP::LibintInterfacer::get().interface(shell);

    BOOST_CHECK(std::equal(ref_libint_shell.alpha.begin(), ref_libint_shell.alpha.end(), libint_shell.alpha.begin(), approx()));  // check Gaussian exponents
    BOOST_CHECK(std::equal(ref_libint_shell.contr[0].coeff.begin(), ref_libint_shell.contr[0].coeff.end(), libint_shell.contr[0].coeff.begin(), approx()));  // check contraction coefficients
    BOOST_CHECK(ref_libint_shell.contr[0].pure == libint_shell.contr[0].pure);  // check if both shells have the same type
    BOOST_CHECK(std::equal(ref_libint_shell.O.begin(), ref_libint_shell.O.end(), libint_shell.O.begin(), approx()));  // check the center
}
