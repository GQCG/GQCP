#define BOOST_TEST_MODULE "LibintCommunicator"


#include "LibintCommunicator.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( atoms_interface ) {

    std::vector<GQCG::Atom> gqcg_atoms = {
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


    /**
     *  Function object for (libint_atom) == (gqcg_atom)
     */
    struct AtomEqual {
        double tolerance;
        explicit AtomEqual(double tolerance) : tolerance (tolerance) {};

        bool operator()(const libint2::Atom& libint_atom, const GQCG::Atom& gqcg_atom) {
            return (gqcg_atom.atomic_number == libint_atom.atomic_number) &&
                   (std::abs(gqcg_atom.x - libint_atom.x) < tolerance) &&
                   (std::abs(gqcg_atom.y - libint_atom.y) < tolerance) &&
                   (std::abs(gqcg_atom.z - libint_atom.z) < tolerance);
        }
    };


    // Check if the interfacing between the Atom types works
    BOOST_CHECK(std::equal(ref_libint_atoms.begin(), ref_libint_atoms.end(), gqcg_atoms.begin(), AtomEqual(1.0e-08)));
}
