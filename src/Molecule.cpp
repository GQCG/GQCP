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
#include "Molecule.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "elements.hpp"
#include "units.hpp"


namespace GQCP {


/*
 *  PRIVATE METHODS
 */

/**
 *  Parse a @param xyz_filename to @return a std::vector<GQCP::Atom>.
 *
 *  The coordinates in the .xyz-file should be in Angstrom: this function converts them immediately to Bohr (a.u.)
 */
std::vector<GQCP::Atom> Molecule::parseXYZFile(const std::string& xyz_filename) {

    // Find the extension of the given path (https://stackoverflow.com/a/51992)
    std::string extension;
    std::string::size_type idx = xyz_filename.rfind('.');

    if (idx != std::string::npos) {
        extension = xyz_filename.substr(idx+1);
    } else {
        throw std::runtime_error("I did not find an extension in your given path.");
    }

    if (!(extension == "xyz")) {
        throw std::runtime_error("You did not provide a .xyz file name");
    }

    // If the xyz_filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream (xyz_filename);
    if (!input_file_stream.good()) {
        throw std::runtime_error("The provided .xyz file name is illegible. Maybe you specified a wrong path?");
    } else {
        // Do the actual parsing
        std::string line;


        // First line is the number of atoms
        std::getline(input_file_stream, line);
        auto number_of_atoms = static_cast<size_t>(std::stoi(line));


        // Second line is empty
        std::getline(input_file_stream, line);


        // Next lines are the atoms
        std::vector<GQCP::Atom> atoms;
        while (std::getline(input_file_stream, line)) {
            std::string symbol;
            double x_angstrom, y_angstrom, z_angstrom;

            std::istringstream iss (line);
            iss >> symbol >> x_angstrom >> y_angstrom >> z_angstrom;

            // Convert the (x,y,z)-coordinates that are in Angstrom to Bohr
            double x_bohr = GQCP::units::angstrom_to_bohr(x_angstrom);
            double y_bohr = GQCP::units::angstrom_to_bohr(y_angstrom);
            double z_bohr = GQCP::units::angstrom_to_bohr(z_angstrom);

            atoms.emplace_back(GQCP::elements::elementToAtomicNumber(symbol), x_bohr, y_bohr, z_bohr);
        }


        if (number_of_atoms > atoms.size()) {
            throw std::invalid_argument("The .xyz-file contains more atoms than specified on its first line.");
        } else {
            return atoms;
        }
    }
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor from a @param atoms: a given std::vector of GQCP::Atoms and a @param molecular_charge
 *      The constructed molecule instance corresponds to an ion:
 *          charge = +1 -> cation (one electron less than the neutral molecule)
 *          charge = 0  -> neutral molecule
 *          charge = -1 -> anion (one electron more than the neutral molecule)
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<GQCP::Atom>& atoms, int molecular_charge) :
    atoms (atoms),
    N (this->calculateTotalNucleicCharge() - molecular_charge)
{
    // Check if the total positive charge is valid, e.g. H^(2+) does not exist
    if (molecular_charge > 0) {
        if (this->calculateTotalNucleicCharge() < molecular_charge) {
            throw std::invalid_argument("You cannot create a molecule with these atoms and this much of a total positive charge.");
        }
    }

    // Check if there are no duplicate atoms
    // The atoms are const, so we will make a copy to check if there are duplicates
    auto atoms_copy = this->atoms;

    // Sort and unique
    std::sort(atoms_copy.begin(), atoms_copy.end());
    auto last_it = std::unique(atoms_copy.begin(), atoms_copy.end());

    // If the iterator returned from unique is the same as the end iterator, there are no duplicate items
    if (last_it != atoms_copy.end()) {
        throw std::invalid_argument("There can't be two equal atoms on the same position.");
    }
}


/**
 *  Constructor from a @param atoms: a given std::vector of GQCP::Atoms
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<GQCP::Atom>& atoms) :
    Molecule (atoms, 0)
{}


/**
 *  Constructor from a given @param xyz_filename and a @param molecular_charge
 *      The constructed molecule instance corresponds to an ion:
 *          charge = +1 -> cation (one electron less than the neutral molecule)
 *          charge = 0  -> neutral molecule
 *          charge = -1 -> anion (one electron more than the neutral molecule)
 *
 *  IMPORTANT!!! The coordinates of the atoms in the .xyz-file should be in Angstrom, but we convert them internally to Bohr
 */
Molecule::Molecule(const std::string& xyz_filename, int molecular_charge) :
    Molecule (Molecule::parseXYZFile(xyz_filename), molecular_charge)
{}


/**
 *  Constructor from a given @param xyz_filename
 *      The constructed molecule instance corresponds to a neutral atom (i.e. N = sum of nucleus charges)
 *
 *  IMPORTANT!!! The coordinates of the atoms in the .xyz-file should be in Angstrom, but we convert them internally to Bohr
 */
Molecule::Molecule(const std::string& xyz_filename) :
    Molecule (xyz_filename, 0)
{}



/*
 *  OPERATORS
 */

/**
 *  @return if this is equal to @param other, within the default GQCP::Atom::tolerance_for_comparison for the coordinates of the atoms
 */
bool Molecule::operator==(const GQCP::Molecule& other) const {

    return this->isEqualTo(other, GQCP::Atom::tolerance_for_comparison);
}


/**
 *  Overloading of operator<< for a GQCP::Molecule to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCP::Molecule& molecule) {

    for (const auto& atom : molecule.atoms) {
        os << atom;
    }

    return os;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return if this is equal to @param other, within the given @param tolerance for the coordinates of the atoms
 */
bool Molecule::isEqualTo(const GQCP::Molecule& other, double tolerance) const {

    if (this->N != other.get_N()) {
        return false;
    }

    // We don't want the order of the atoms to matter in a GQCP::Molecule comparison
    // We have implemented a custom GQCP::Atom::operator< so we can sort std::vectors of GQCP::Atoms
    // Make a copy of the atoms because std::sort modifies
    auto this_atoms = this->atoms;
    auto other_atoms = other.atoms;


    // Make lambda expressions for the comparators, since we want to hand over the tolerance argument
    auto smaller_than_atom = [this, tolerance](const GQCP::Atom& lhs, const GQCP::Atom& rhs) { return lhs.isSmallerThan(rhs, tolerance); };
    auto equal_atom = [this, tolerance](const GQCP::Atom& lhs, const GQCP::Atom& rhs) { return lhs.isEqualTo(rhs, tolerance); };


    std::sort(this_atoms.begin(), this_atoms.end(), smaller_than_atom);
    std::sort(other_atoms.begin(), other_atoms.end(), smaller_than_atom);

    return std::equal(this_atoms.begin(), this_atoms.end(), other_atoms.begin(), equal_atom);
}


/**
 *  @return the sum of all the charges of the nuclei
 */
size_t Molecule::calculateTotalNucleicCharge() const {

    size_t nucleic_charge = 0;

    for (const auto& atom : this->atoms) {
        nucleic_charge += atom.atomic_number;
    }

    return nucleic_charge;
}


/**
 * @return the distance between two the two atoms at @param index1 and @param index2
 */
double Molecule::calculateInternuclearDistance(size_t index1, size_t index2) const {

    // Check if the indices are within bounds
    if (index1 > this->atoms.size() || index2 > this->atoms.size()) {
        throw std::invalid_argument("At least one of the given indices is out of bounds.");
    }

    return this->atoms[index1].calculateDistance(this->atoms[index2]);
}


/**
 *  @return the internuclear repulsion energy due to the nuclear framework
 */
double Molecule::calculateInternuclearRepulsionEnergy() const {

    double internuclear_repulsion_energy = 0.0;

    // Sum over every unique nucleus pair
    auto natoms = this->numberOfAtoms();
    for (size_t i = 0; i < natoms; i++) {
        for (size_t j = i + 1; j < natoms; j++ ) {
            const auto atom1 = this->atoms[i];
            const auto atom2 = this->atoms[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
            internuclear_repulsion_energy += atom1.atomic_number * atom2.atomic_number / this->calculateInternuclearDistance(i, j);
        }
    }

    return internuclear_repulsion_energy;
}


}  // namespace GQCP
