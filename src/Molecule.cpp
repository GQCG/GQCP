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
#include "Molecule.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "elements.hpp"
#include "units.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */

/**
 *  @param atoms        the atoms that make up the molecule, with coordinates in bohr
 *  @param charge       the charge of the molecule:
 *                          +1 -> cation (one electron less than the neutral molecule)
 *                           0 -> neutral molecule
 *                          -1 -> anion (one electron more than the neutral molecule)
 */
Molecule::Molecule(const std::vector<Atom>& atoms, int charge) :
    atoms (atoms),
    N (this->calculateTotalNucleicCharge() - charge)
{
    // Check if the total positive charge is valid, e.g. H^(2+) does not exist
    if (charge > 0) {
        if (this->calculateTotalNucleicCharge() < charge) {
            throw std::invalid_argument("You cannot create a molecule with these atoms and this much of a total positive charge.");
        }
    }

    // Check if there are no duplicate atoms
    std::vector<Atom> atoms_copy = this->atoms;

    // Sort and unique
    std::sort(atoms_copy.begin(), atoms_copy.end());
    auto last_it = std::unique(atoms_copy.begin(), atoms_copy.end());

    // If the iterator returned from unique is the same as the end iterator, there are no duplicate items
    if (last_it != atoms_copy.end()) {
        throw std::invalid_argument("There can't be two equal atoms on the same position.");
    }
}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
 *
 *  @return a vector of Atoms that are in the given xyz-file
 */
Molecule Molecule::Readxyz(const std::string& xyz_filename, int charge) {

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
        std::vector<Atom> atoms;
        atoms.reserve(number_of_atoms);

        while (std::getline(input_file_stream, line)) {
            std::string symbol;
            double x_angstrom, y_angstrom, z_angstrom;

            std::istringstream iss (line);
            iss >> symbol >> x_angstrom >> y_angstrom >> z_angstrom;

            // Convert the (x,y,z)-coordinates that are in Angstrom to Bohr
            double x_bohr = units::angstrom_to_bohr(x_angstrom);
            double y_bohr = units::angstrom_to_bohr(y_angstrom);
            double z_bohr = units::angstrom_to_bohr(z_angstrom);

            atoms.emplace_back(elements::elementToAtomicNumber(symbol), x_bohr, y_bohr, z_bohr);
        }


        if (number_of_atoms > atoms.size()) {
            throw std::invalid_argument("The .xyz-file contains more atoms than specified on its first line.");
        } else {
            return Molecule(atoms, charge);
        }
    }
}


/**
 *  @param n            the number of H atoms
 *  @param spacing      the internuclear spacing in bohr
 *  @param charge       the total charge
 *
 *  @return a charged H-chain with equal internuclear spacing
 */
Molecule Molecule::HChain(size_t n, double spacing, int charge) {

    if (n == 0) {
        throw std::invalid_argument("Can not create a H-chain consisting of zero atoms.");
    }

    if (spacing < 0.0) {
        throw std::invalid_argument("Can't have a negative spacing.");
    }


    std::vector<Atom> h_chain;

    // Put all H-atoms on a line on the x-axis: the first H is on the origin
    double x = 0.0;  // the current x-coordinate
    for (size_t i = 0; i < n; i++) {
        h_chain.emplace_back(1, x, 0.0, 0.0);

        x += spacing;  // proceed to the next H-atom
    }

    return Molecule(h_chain, charge);
}


/**
 *  @param n        the number of H2-molecules
 *  @param a        the internuclear distance in bohr
 *  @param b        the intermolecular distance in bohr
 *  @param charge   the total charge
 *
 *  @return a charged H2-chain
 */
Molecule Molecule::H2Chain(size_t n, double a, double b, int charge) {

    if (n == 0) {
        throw std::invalid_argument("Can not create a H2-chain consisting of zero H2-molecules.");
    }

    if ((a < 0.0) || (b < 0.0)) {
        throw std::invalid_argument("Can't have a negative spacing.");
    }


    std::vector<Atom> h_chain;

    // Put all H-atoms on a line on the x-axis: the first H is on the origin
    double x = 0.0;  // the current x-coordinate
    for (size_t i = 0; i < n; i++) {

        h_chain.emplace_back(1, x, 0.0, 0.0);  // the first H-atom
        x += a;  // add internuclear distance (we're within a H2-molecule)
        h_chain.emplace_back(1, x, 0.0, 0.0);  // the second H-atom

        x += b;  // proceed to the next H2-molecule
    }

    return Molecule(h_chain, charge);
}



/*
 *  OPERATORS
 */

/**
 *  @param other        the other molecule
 *
 *  @return if this molecule is equal to the other, within the default Atom::tolerance_for_comparison for the coordinates of the atoms
 */
bool Molecule::operator==(const Molecule& other) const {

    return this->isEqualTo(other, Atom::tolerance_for_comparison);
}


/**
 *  @param os           the output stream which the molecule should be concatenated to
 *  @param molecule     the molecule that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {

    for (const auto& atom : molecule.atoms) {
        os << atom;
    }

    return os;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other molecule
 *  @param tolerance    the tolerance for the coordinates of the atoms
 *
 *  @return if this is equal to the other, within the given tolerance
 */
bool Molecule::isEqualTo(const Molecule& other, double tolerance) const {

    if (this->N != other.get_N()) {
        return false;
    }

    // We don't want the order of the atoms to matter in a Molecule comparison
    // We have implemented a custom Atom::operator< so we can sort std::vectors of Atoms
    // Make a copy of the atoms because std::sort modifies
    auto this_atoms = this->atoms;
    auto other_atoms = other.atoms;


    // Make lambda expressions for the comparators, since we want to hand over the tolerance argument
    auto smaller_than_atom = [this, tolerance](const Atom& lhs, const Atom& rhs) { return lhs.isSmallerThan(rhs, tolerance); };
    auto equal_atom = [this, tolerance](const Atom& lhs, const Atom& rhs) { return lhs.isEqualTo(rhs, tolerance); };


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
 *  @param index1   the index of the first atom
 *  @param index2   the index of the second atom
 *
 *  @return the distance between the two atoms at index1 and index2 in bohr
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


/**
 *  @return the electrical dipole moment generated by the nuclear framework
 */
Eigen::Vector3d Molecule::calculateNuclearDipoleMoment() const {

    Eigen::Vector3d m = Eigen::Vector3d::Zero();

    for (const auto& atom : this->atoms) {
        m += atom.atomic_number * atom.position;
    }

    return m;
}


}  // namespace GQCP
