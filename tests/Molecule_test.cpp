#define BOOST_TEST_MODULE "Molecule"


#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( constructor_atoms_charge ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCG::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    // Check if we can create any anion
    GQCG::Molecule molecule2 (atoms, -2);


    // Check if we can't create a cation with charge larger than the nucleic charge
    BOOST_CHECK_NO_THROW(GQCG::Molecule (atoms, +3));
    BOOST_CHECK_THROW(GQCG::Molecule (atoms, +11), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_atoms ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCG::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    GQCG::Molecule molecule (atoms);
}


BOOST_AUTO_TEST_CASE ( calculateTotalNucleicCharge ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCG::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    GQCG::Molecule molecule (atoms);
    BOOST_CHECK_EQUAL(molecule.calculateTotalNucleicCharge(), 10);
}


BOOST_AUTO_TEST_CASE ( parseXYZFile ) {

    // Make sure we get an error when a nonsense path is given for the .xyz file name
    BOOST_REQUIRE_THROW(GQCG::Molecule ("this is a nonsense data path"), std::runtime_error);

    // Make sure we get an error when a path with a wrong extension is given
    BOOST_REQUIRE_THROW(GQCG::Molecule ("../tests/ref_data/nuclear.data"), std::runtime_error);

    // Make sure we don't get an error when a correct path is given
    BOOST_REQUIRE_NO_THROW(GQCG::Molecule ("../tests/data/h2o.xyz"));
}


BOOST_AUTO_TEST_CASE ( molecule_ion_constructor ) {

    // Create some Molecule objects
    const std::string xyzfilename = "../tests/data/h2o.xyz";
    GQCG::Molecule water (xyzfilename);
    GQCG::Molecule water_anion (xyzfilename, -1);
    GQCG::Molecule water_neutral (xyzfilename, 0);
    GQCG::Molecule water_cation (xyzfilename, +1);

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.get_N(), 10);
    BOOST_CHECK_EQUAL(water_anion.get_N(), 11);
    BOOST_CHECK_EQUAL(water_neutral.get_N(), 10);
    BOOST_CHECK_EQUAL(water_cation.get_N(), 9);
}


BOOST_AUTO_TEST_CASE ( xyz_filename_constructor ) {



}


