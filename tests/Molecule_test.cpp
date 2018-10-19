// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "Molecule"


#include "Molecule.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise clang++ will complain


BOOST_AUTO_TEST_CASE ( constructor_atoms_charge ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCP::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    // Check if we can create any anion
    GQCP::Molecule molecule2 (atoms, -2);


    // Check if we can't create a cation with charge larger than the nucleic charge
    BOOST_CHECK_NO_THROW(GQCP::Molecule (atoms, +3));
    BOOST_CHECK_THROW(GQCP::Molecule (atoms, +11), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( constructor_atoms ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCP::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    GQCP::Molecule molecule (atoms);
}


BOOST_AUTO_TEST_CASE ( duplicate_atoms_constructor ) {

    // Make some atoms
    GQCP::Atom atom1 {1, 0.0, 0.0, 0.0};
    GQCP::Atom atom2 {1, 1.0, 0.0, 0.0};

    std::vector<GQCP::Atom> atoms1 {atom1, atom1};
    std::vector<GQCP::Atom> atoms2 {atom1, atom2};


    // Check if we can't create a Molecule with duplicate atoms
    BOOST_CHECK_THROW(GQCP::Molecule molecule (atoms1), std::invalid_argument);

    // Check if a correct argument doesn't throw
    BOOST_CHECK_NO_THROW(GQCP::Molecule molecule (atoms2));
}


BOOST_AUTO_TEST_CASE ( calculateTotalNucleicCharge ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCP::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };

    GQCP::Molecule molecule (atoms);
    BOOST_CHECK_EQUAL(molecule.calculateTotalNucleicCharge(), 10);
}


BOOST_AUTO_TEST_CASE ( parseXYZFile ) {

    // Make sure we get an error when a nonsense path is given for the .xyz file name
    BOOST_REQUIRE_THROW(GQCP::Molecule ("this is a nonsense data path"), std::runtime_error);

    // Make sure we get an error when a path with a wrong extension is given
    BOOST_REQUIRE_THROW(GQCP::Molecule ("../tests/ref_data/nuclear.data"), std::runtime_error);

    // Make sure we don't get an error when a correct path is given
    BOOST_REQUIRE_NO_THROW(GQCP::Molecule ("../tests/data/h2o.xyz"));
}


BOOST_AUTO_TEST_CASE ( molecule_ion_constructor ) {

    // Create some Molecule objects
    const std::string xyzfilename = "../tests/data/h2o.xyz";
    GQCP::Molecule water (xyzfilename);
    GQCP::Molecule water_anion (xyzfilename, -1);
    GQCP::Molecule water_neutral (xyzfilename, 0);
    GQCP::Molecule water_cation (xyzfilename, +1);

    // Test the number of electrons created by the constructor
    BOOST_CHECK_EQUAL(water.get_N(), 10);
    BOOST_CHECK_EQUAL(water_anion.get_N(), 11);
    BOOST_CHECK_EQUAL(water_neutral.get_N(), 10);
    BOOST_CHECK_EQUAL(water_cation.get_N(), 9);
}


BOOST_AUTO_TEST_CASE ( Molecule_operator_ostream ) {

    std::vector<GQCP::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };
    GQCP::Molecule molecule (atoms);

    std::cout << molecule << std::endl;
}


BOOST_AUTO_TEST_CASE ( Molecule_isEqualTo ) {

    // Create some Atoms and Molecules
    GQCP::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCP::Atom atom2 {2, 0.0, 0.1, 0.2};
    GQCP::Atom atom3 {3, 0.0, 0.1, 0.2};
    GQCP::Atom atom4 {4, 0.1, 0.2, 0.3};
    GQCP::Atom atom5 {3, 0.1, 0.2, 0.3};

    GQCP::Molecule molecule1 {{atom1, atom2, atom3}};
    GQCP::Molecule molecule2 {{atom1, atom2, atom3}};
    GQCP::Molecule molecule3 {{atom1, atom2, atom3}, -1};
    GQCP::Molecule molecule4 {{atom1, atom2, atom5}};
    GQCP::Molecule molecule5 {{atom1, atom3, atom2}};
    GQCP::Molecule molecule6 {{atom1, atom2, atom3, atom4}};

    // Check if they're equal
    BOOST_CHECK(molecule1.isEqualTo(molecule2));

    // Check if a different charge but same atoms causes inequality
    BOOST_CHECK(!(molecule1.isEqualTo(molecule3)));

    // Check if different atoms but an equal total charge cause inequality
    BOOST_CHECK(!(molecule1.isEqualTo(molecule4)));

    // Check if a different ordering doesn't cause inequality
    BOOST_CHECK(molecule1.isEqualTo(molecule5));

    // Check if a different number of atoms causes inequality
    BOOST_CHECK(!(molecule1.isEqualTo(molecule6)));


    // Check if the tolerance argument works
    BOOST_CHECK(molecule1.isEqualTo(molecule4, 0.2));
}


BOOST_AUTO_TEST_CASE ( Molecule_operator_equals ) {

    // Create some Atoms and Molecules
    GQCP::Atom atom1 {1, 0.0, 0.1, 0.2};
    GQCP::Atom atom2 {2, 0.0, 0.1, 0.2};
    GQCP::Atom atom3 {3, 0.0, 0.1, 0.2};
    GQCP::Atom atom4 {3, 0.1, 0.2, 0.3};

    GQCP::Molecule molecule1 {{atom1, atom2, atom3}};
    GQCP::Molecule molecule2 {{atom1, atom2, atom3}};
    GQCP::Molecule molecule3 {{atom1, atom2, atom4}};


    // Check if we can call operator==
    BOOST_CHECK(molecule1 == molecule2);
    BOOST_CHECK(!(molecule2 == molecule3));
}


BOOST_AUTO_TEST_CASE ( xyz_filename_constructor ) {

    std::vector<GQCP::Atom> atoms {
        {8,  0.0,     -0.143222, 0.0},
        {1,  1.63803,  1.13656,  0.0},
        {1, -1.63803,  1.13656,  0.0}
    };
    GQCP::Molecule molecule_atoms (atoms);

    GQCP::Molecule molecule_xyz ("../tests/data/h2o.xyz");

    // Check if the conversion from Bohr to Angstrom is correct
    BOOST_CHECK(molecule_atoms.isEqualTo(molecule_xyz, 1.0e-05));
}


BOOST_AUTO_TEST_CASE ( calculateInternuclearDistance ) {

    // Create a fictitious molecule from some Atoms (charge, x, y ,z)
    std::vector<GQCP::Atom> atoms = {
        {1, 0, 3, 0},
        {2, 0, 0, 4},
        {3, 3, 0, 0},
        {4, 0, 0, 5}
    };
    GQCP::Molecule molecule (atoms);


    // Check if we get throws when the indices are out of bounds
    BOOST_CHECK_THROW(molecule.calculateInternuclearDistance(0, 5), std::invalid_argument);
    BOOST_CHECK_THROW(molecule.calculateInternuclearDistance(8, 2), std::invalid_argument);

    // Check if we don't get throws when the indices behave correctly
    BOOST_CHECK_NO_THROW(molecule.calculateInternuclearDistance(0, 0));

    // Check if the function works
    BOOST_CHECK(std::abs(molecule.calculateInternuclearDistance(1, 3) - 1) < 1.0e-12);
}


BOOST_AUTO_TEST_CASE ( methods_h2 ) {

    // We have reference internuclear repulsion energy from HORTON
    double ref_internuclear_repulsion_energy = 0.714285658963;

    // Create the hydrogen gas molecule
    GQCP::Molecule h2 ("../tests/data/h2_szabo.xyz");

    // Test the basic methods
    BOOST_CHECK_EQUAL(h2.numberOfAtoms(), 2);
    BOOST_CHECK_EQUAL(h2.calculateTotalNucleicCharge(), 2);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(h2.calculateInternuclearRepulsionEnergy() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}


BOOST_AUTO_TEST_CASE ( methods_water ) {

    // We have reference internuclear repulsion energy from HORTON
    double ref_internuclear_repulsion_energy = 8.00236693455;

    // Create the water molecule
    GQCP::Molecule water ("../tests/data/h2o.xyz");

    // Test the basic methods
    BOOST_CHECK_EQUAL(water.numberOfAtoms(), 3);
    BOOST_CHECK_EQUAL(water.calculateTotalNucleicCharge(), 10);

    // Test the calculation of the nuclear repulsion energy
    BOOST_CHECK(std::abs(water.calculateInternuclearRepulsionEnergy() - ref_internuclear_repulsion_energy) < 1.0e-07);  // reference data from horton
}
