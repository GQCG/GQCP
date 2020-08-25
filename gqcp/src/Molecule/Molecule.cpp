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

#include "Molecule/Molecule.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param nuclear_framework        the nuclear framework that makes up the molecule, with coordinates in bohr
 *  @param charge                   the charge of the molecule:
 *                                      +1 -> cation (one electron less than the neutral molecule)
 *                                       0 -> neutral molecule
 *                                      -1 -> anion (one electron more than the neutral molecule)
 */
Molecule::Molecule(const NuclearFramework& nuclear_framework, const int charge) :
    nuclear_framework {nuclear_framework},
    N {nuclear_framework.totalNucleicCharge() - charge} {

    // Check if the total positive charge is valid, e.g. H^(2+) does not exist
    if (charge > 0) {
        if (this->totalNucleicCharge() < charge) {  // OK to compare size_t and int like this
            throw std::invalid_argument("Molecule::Molecule(NuclearFramework,int): You cannot create a molecule with these nuclei and this much of a total positive charge.");
        }
    }
}


/**
 *  @param nuclei       the nuclei that make up the molecule, with coordinates in bohr
 *  @param charge       the charge of the molecule:
 *                          +1 -> cation (one electron less than the neutral molecule)
 *                           0 -> neutral molecule
 *                          -1 -> anion (one electron more than the neutral molecule)
 */
Molecule::Molecule(const std::vector<Nucleus>& nuclei, const int charge) :
    Molecule(NuclearFramework(nuclei), charge) {}


/*
 *  NAMED CONSTRUCTORS
 */

/**
 *  @param n            the number of H nuclei
 *  @param spacing      the internuclear spacing in bohr
 *  @param charge       the total charge
 *
 *  @return a charged H-chain with equal internuclear spacing
 */
Molecule Molecule::HChain(const size_t n, const double spacing, const int charge, const CartesianDirection axis) {

    return Molecule(NuclearFramework::HChain(n, spacing, axis), charge);
}


/**
 *  @param n        the number of H2-molecules
 *  @param a        the internuclear distance in bohr
 *  @param b        the intermolecular distance in bohr
 *  @param charge   the total charge
 *
 *  @return a charged H2-chain
 */
Molecule Molecule::H2Chain(const size_t n, const double a, const double b, const int charge, const CartesianDirection axis) {

    return Molecule(NuclearFramework::H2Chain(n, a, b, axis), charge);
}


/**
 *  @param n                the number of hydrogens
 *  @param distance         the distance (in bohr) between neighbouring hydrogen atoms
 * 
 *  @return a regular H-ring where neighbouring hydrogens are separated by the given distance
 */
Molecule Molecule::HRingFromDistance(const size_t n, const double distance, const int charge) {

    return Molecule(NuclearFramework::HRingFromDistance(n, distance), charge);
}


/**
 *  @param n                the number of hydrogens
 *  @param radius           the radius (in bohr) of the circumscribed circle
 * 
 *  @return a regular H-ring whose hydrogens are on the circle with the given radius
 */
Molecule Molecule::HRingFromRadius(const size_t n, const double radius, const int charge) {

    return Molecule(NuclearFramework::HRingFromRadius(n, radius), charge);
}


/**
 *  Construct a molecule based on the content of a given .xyz-file. In an .xyz-file, the molecular coordinates are in Angstrom
 *
 *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
 *  @param charge       the charge of the molecule:
 *                          +1 -> cation (one electron less than the neutral molecule)
 *                           0 -> neutral molecule
 *                          -1 -> anion (one electron more than the neutral molecule)
 */
Molecule Molecule::ReadXYZ(const std::string& xyz_filename, const int charge) {

    return Molecule(NuclearFramework::ReadXYZ(xyz_filename), charge);
}


/*
 *  OPERATORS
 */

/**
 *  @param os           the output stream which the molecule should be concatenated to
 *  @param molecule     the molecule that should be concatenated to the output stream
 *
 *  @return the updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {

    os << molecule.description();

    return os;
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the charge of this molecule (in a.u.)
 */
int Molecule::charge() const {
    return this->totalNucleicCharge() - this->numberOfElectrons();
}


}  // namespace GQCP
