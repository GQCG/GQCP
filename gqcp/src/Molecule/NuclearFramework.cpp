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
#include "Molecule/NuclearFramework.hpp"

#include "Molecule/elements.hpp"
#include "Utilities/units.hpp"
#include "Utilities/miscellaneous.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param nuclei           the nuclei of the nuclear framework
 */
NuclearFramework::NuclearFramework(const std::vector<Nucleus>& nuclei) :
    nuclei (nuclei)
{
    // Check if there are no duplicate nuclei
    std::vector<Nucleus> nuclei_copy = this->nuclei;

    // Sort and unique
    std::sort(nuclei_copy.begin(), nuclei_copy.end(), Nucleus::sortComparer());
    auto last_it = std::unique(nuclei_copy.begin(), nuclei_copy.end(), Nucleus::equalityComparer());

    // If the iterator returned from unique is the same as the end iterator, there are no duplicate items
    if (last_it != nuclei_copy.end()) {
        throw std::invalid_argument("NuclearFramework::NuclearFramework(std::vector<Nucleus>): There can't be two nuclei on the same position.");
    }
}



/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  Construct a nuclear framework based on the content of a given .xyz-file. In an .xyz-file, the nuclear coordinates are in Angstrom
 *
 *  @param xyz_filename     the .xyz-file that contains the nuclear coordinates in Angstrom
 */
NuclearFramework NuclearFramework::ReadXYZ(const std::string& xyz_filename) {

    // Check the file name extension and open the file
    std::ifstream input_file_stream = validateAndOpen(xyz_filename, "xyz");


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
        throw std::invalid_argument("NuclearFramework::ReadXYZ(const std::string&): The .xyz-file contains more nuclei than specified on its first line.");
    } else {
        return NuclearFramework(nuclei);
    }
}


/**
 *  @param n            the number of H nuclei
 *  @param spacing      the internuclear spacing in bohr
 *
 *  @return a H-chain with equal internuclear spacing
 */
NuclearFramework NuclearFramework::HChain(const size_t n, const double spacing, CartesianDirection axis) {

    if (n == 0) {
        throw std::invalid_argument("NuclearFramework::HChain(const size_t, const double): Can not create a H-chain consisting of zero nuclei.");
    }

    if (spacing < 0.0) {
        throw std::invalid_argument("NuclearFramework::HChain(const size_t, const double): Can't have a negative spacing.");
    }


    std::vector<Nucleus> h_chain;

    // Put all H-nuclei on a line on the given axis: the first H is on the origin
    double coordinate {0.0};  // the current x-coordinate
    Vector<double, 3> position = Vector<double, 3>::Zero(3);
    for (size_t i = 0; i < n; i++) {
        position(axis) = coordinate;  // put a value on index 0, 1, or 2

        h_chain.emplace_back(1, position);

        coordinate += spacing;  // proceed to the next H-atom
    }

    return NuclearFramework(h_chain);
}


/**
 *  @param n        the number of H2-molecules
 *  @param a        the internuclear distance in bohr
 *  @param b        the intermolecular distance in bohr
 *
 *  @return a H2-chain with the specified internuclear and intermolecular distances
 */
NuclearFramework NuclearFramework::H2Chain(const size_t n, const double a, const double b, CartesianDirection axis) {

    if (n == 0) {
        throw std::invalid_argument("NuclearFramework::H2Chain(const size_t, const double, const double): Can not create a H2-chain consisting of zero H2-molecules.");
    }

    if ((a < 0.0) || (b < 0.0)) {
        throw std::invalid_argument("NuclearFramework::H2Chain(const size_t, const double, const double): Can't have a negative spacing.");
    }


    std::vector<Nucleus> h_chain;

    // Put all H-nuclei on a line on the given: the first H is on the origin
    double coordinate {0.0};  // the current coordinate
    Vector<double, 3> position = Vector<double, 3>::Zero(3);

    for (size_t i = 0; i < n; i++) {

        position(axis) = coordinate;  // put a value on index 0, 1, or 2
        h_chain.emplace_back(1, position);  // the first H-atom

        coordinate += a;  // add internuclear distance (we're within a H2-molecule)

        position(axis) = coordinate;  // put a value on index 0, 1, or 2
        h_chain.emplace_back(1, position);  // the second H-atom

        coordinate += b;  // proceed to the next H2-molecule
    }

    return NuclearFramework(h_chain);
}


/*
 *  OPERATORS
 */

/**
 *  @param os                       the output stream which the nuclear framework should be concatenated to
 *  @param nuclear_framework        the nuclear framework that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const NuclearFramework& nuclear_framework) {

    for (const auto& nucleus : nuclear_framework.nucleiAsVector()) {
        os << nucleus;
    }

    return os;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the sum of all the charges of the nuclei
 */
size_t NuclearFramework::totalNucleicCharge() const {

    size_t value {0};

    for (const auto& nucleus : this->nuclei) {
        value += nucleus.charge();
    }

    return value;
}


/**
 *  @param index1   the index of the first nucleus
 *  @param index2   the index of the second nucleus
 *
 *  @return the distance between the two nuclei at index1 and index2 in bohr
 */
double NuclearFramework::internuclearDistance(size_t index1, size_t index2) const {

    // Check if the indices are within bounds
    if (index1 > this->nuclei.size() || index2 > this->nuclei.size()) {
        throw std::invalid_argument("NuclearFramework::internuclearDistance(const size_t, const size_t): At least one of the given indices is out of bounds.");
    }

    return this->nuclei[index1].calculateDistance(this->nuclei[index2]);
}



}  // namespace GQCP
