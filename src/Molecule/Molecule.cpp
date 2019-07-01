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
#include "Molecule/Molecule.hpp"

#include "elements.hpp"
#include "units.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>



namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param nuclei       the nuclei that make up the molecule, with coordinates in bohr
 *  @param charge       the charge of the molecule:
 *                          +1 -> cation (one electron less than the neutral molecule)
 *                           0 -> neutral molecule
 *                          -1 -> anion (one electron more than the neutral molecule)
 */
Molecule::Molecule(const std::vector<Nucleus>& nuclei, int charge) :
    nuclei (nuclei),
    N (this->calculateTotalNucleicCharge() - charge)
{
    // Check if the total positive charge is valid, e.g. H^(2+) does not exist
    if (charge > 0) {
        if (this->calculateTotalNucleicCharge() < charge) {
            throw std::invalid_argument("Molecule::Molecule(): You cannot create a molecule with these nuclei and this much of a total positive charge.");
        }
    }

    // Check if there are no duplicate nuclei
    std::vector<Nucleus> nuclei_copy = this->nuclei;

    // Sort and unique
    std::sort(nuclei_copy.begin(), nuclei_copy.end());
    auto last_it = std::unique(nuclei_copy.begin(), nuclei_copy.end());

    // If the iterator returned from unique is the same as the end iterator, there are no duplicate items
    if (last_it != nuclei_copy.end()) {
        throw std::invalid_argument("Molecule::Molecule(): There can't be two equal nuclei on the same position.");
    }
}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Construct a molecule based on the content of a given .xyz-file. In an .xyz-file, the molecular coordinates are in Angstrom
 *
 *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
 *  @param charge       the charge of the molecule:
 *                          +1 -> cation (one electron less than the neutral molecule)
 *                           0 -> neutral molecule
 *                          -1 -> anion (one electron more than the neutral molecule)
 */
Molecule Molecule::Readxyz(const std::string& xyz_filename, int charge) {

    // Find the extension of the given path (https://stackoverflow.com/a/51992)
    std::string extension;
    std::string::size_type idx = xyz_filename.rfind('.');

    if (idx != std::string::npos) {
        extension = xyz_filename.substr(idx+1);
    } else {
        throw std::invalid_argument("Molecule::Readxyz(const std::string&, int): I did not find an extension in your given path.");
    }

    if (!(extension == "xyz")) {
        throw std::invalid_argument("Molecule::Readxyz(const std::string&, int): You did not provide a .xyz file name");
    }

    // If the xyz_filename isn't properly converted into an input file stream, we assume the user supplied a wrong file
    std::ifstream input_file_stream (xyz_filename);
    if (!input_file_stream.good()) {
        throw std::invalid_argument("Molecule::Readxyz(const std::string&, int): The provided .xyz file name is illegible. Maybe you specified a wrong path?");
    } else {
        // Do the actual parsing
        std::string line;


        // First line is the number of nuclei
        std::getline(input_file_stream, line);
        auto number_of_nuclei = static_cast<size_t>(std::stoi(line));


        // Second line is empty
        std::getline(input_file_stream, line);


        // Next lines are the nuclei
        std::vector<Nucleus> nuclei;
        nuclei.reserve(number_of_nuclei);

        while (std::getline(input_file_stream, line)) {
            std::string symbol;
            double x_angstrom, y_angstrom, z_angstrom;

            std::istringstream iss (line);
            iss >> symbol >> x_angstrom >> y_angstrom >> z_angstrom;

            // Convert the (x,y,z)-coordinates that are in Angstrom to Bohr
            double x_bohr = units::angstrom_to_bohr(x_angstrom);
            double y_bohr = units::angstrom_to_bohr(y_angstrom);
            double z_bohr = units::angstrom_to_bohr(z_angstrom);

            nuclei.emplace_back(elements::elementToAtomicNumber(symbol), x_bohr, y_bohr, z_bohr);
        }


        if (number_of_nuclei > nuclei.size()) {
            throw std::invalid_argument("Molecule::Readxyz(const std::string&, int): The .xyz-file contains more nuclei than specified on its first line.");
        } else {
            return Molecule(nuclei, charge);
        }
    }
}


/**
 *  @param n            the number of H nuclei
 *  @param spacing      the internuclear spacing in bohr
 *  @param charge       the total charge
 *
 *  @return a charged H-chain with equal internuclear spacing
 */
Molecule Molecule::HChain(size_t n, double spacing, int charge) {

    if (n == 0) {
        throw std::invalid_argument("Molecule::HChain(size_t, double, int): Can not create a H-chain consisting of zero nuclei.");
    }

    if (spacing < 0.0) {
        throw std::invalid_argument("Molecule::HChain(size_t, double, int): Can't have a negative spacing.");
    }


    std::vector<Nucleus> h_chain;

    // Put all H-nuclei on a line on the x-axis: the first H is on the origin
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
        throw std::invalid_argument("Molecule::H2Chain(size_t, double, double, int): Can not create a H2-chain consisting of zero H2-molecules.");
    }

    if ((a < 0.0) || (b < 0.0)) {
        throw std::invalid_argument("Molecule::H2Chain(size_t, double, double, int): Can't have a negative spacing.");
    }


    std::vector<Nucleus> h_chain;

    // Put all H-nuclei on a line on the x-axis: the first H is on the origin
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
 *  @return if this molecule is equal to the other, within the default Nucleus::tolerance_for_comparison for the coordinates of the nuclei
 */
bool Molecule::operator==(const Molecule& other) const {

    return this->isEqualTo(other, Nucleus::tolerance_for_comparison);
}


/**
 *  @param os           the output stream which the molecule should be concatenated to
 *  @param molecule     the molecule that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {

    for (const auto& nucleus : molecule.nuclei) {
        os << nucleus;
    }

    return os;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param other        the other molecule
 *  @param tolerance    the tolerance for the coordinates of the nuclei
 *
 *  @return if this is equal to the other, within the given tolerance
 */
bool Molecule::isEqualTo(const Molecule& other, double tolerance) const {

    if (this->N != other.get_N()) {
        return false;
    }

    // We don't want the order of the nuclei to matter in a Molecule comparison
    // We have implemented a custom Nucleus::operator< so we can sort std::vectors of Nuclei
    // Make a copy of the nuclei because std::sort modifies
    auto this_nuclei = this->nuclei;
    auto other_nuclei = other.nuclei;


    // Make lambda expressions for the comparators, since we want to hand over the tolerance argument
    auto smaller_than_nucleus = [this, tolerance](const Nucleus& lhs, const Nucleus& rhs) { return lhs.isSmallerThan(rhs, tolerance); };
    auto equal_nucleus = [this, tolerance](const Nucleus& lhs, const Nucleus& rhs) { return lhs.isEqualTo(rhs, tolerance); };


    std::sort(this_nuclei.begin(), this_nuclei.end(), smaller_than_nucleus);
    std::sort(other_nuclei.begin(), other_nuclei.end(), smaller_than_nucleus);

    return std::equal(this_nuclei.begin(), this_nuclei.end(), other_nuclei.begin(), equal_nucleus);
}


/**
 *  @return the sum of all the charges of the nuclei
 */
size_t Molecule::calculateTotalNucleicCharge() const {

    size_t nucleic_charge = 0;

    for (const auto& nucleus : this->nuclei) {
        nucleic_charge += nucleus.atomic_number;
    }

    return nucleic_charge;
}


/**
 *  @param index1   the index of the first nucleus
 *  @param index2   the index of the second nucleus
 *
 *  @return the distance between the two nuclei at index1 and index2 in bohr
 */
double Molecule::calculateInternuclearDistance(size_t index1, size_t index2) const {

    // Check if the indices are within bounds
    if (index1 > this->nuclei.size() || index2 > this->nuclei.size()) {
        throw std::invalid_argument("Molecule::calculateInternuclearDistance(size_t, size_t): At least one of the given indices is out of bounds.");
    }

    return this->nuclei[index1].calculateDistance(this->nuclei[index2]);
}


/**
 *  @return the internuclear repulsion energy due to the nuclear framework
 */
double Molecule::calculateInternuclearRepulsionEnergy() const {

    double internuclear_repulsion_energy = 0.0;

    // Sum over every unique nucleus pair
    auto n_nuclei = this->numberOfAtoms();
    for (size_t i = 0; i < n_nuclei; i++) {
        for (size_t j = i + 1; j < n_nuclei; j++ ) {
            const auto nucleus1 = this->nuclei[i];
            const auto nucleus2 = this->nuclei[j];

            // The internuclear repulsion energy (Coulomb) for every nucleus pair is Z1 * Z2 / |R1 - R2|
            internuclear_repulsion_energy += nucleus1.atomic_number * nucleus2.atomic_number / this->calculateInternuclearDistance(i, j);
        }
    }

    return internuclear_repulsion_energy;
}


/**
 *  @return the electrical dipole moment generated by the nuclear framework
 */
Vector<double, 3> Molecule::calculateNuclearDipoleMoment() const {

    Vector<double, 3> m = Vector<double, 3>::Zero();

    for (const auto& nucleus : this->nuclei) {
        m += nucleus.atomic_number * nucleus.position;
    }

    return m;
}


}  // namespace GQCP
