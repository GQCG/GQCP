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
    nuclear_framework (nuclear_framework),
    N (nuclear_framework.totalNucleicCharge() - charge)
{
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
Molecule::Molecule(const std::vector<Nucleus>& nuclei, int charge) :
    Molecule(NuclearFramework(nuclei), charge)
{}


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
Molecule Molecule::ReadXYZ(const std::string& xyz_filename, int charge) {

    return Molecule(NuclearFramework::ReadXYZ(xyz_filename), charge);
}


/**
 *  @param n            the number of H nuclei
 *  @param spacing      the internuclear spacing in bohr
 *  @param charge       the total charge
 *
 *  @return a charged H-chain with equal internuclear spacing
 */
Molecule Molecule::HChain(size_t n, double spacing, int charge, CartesianDirection axis) {

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
Molecule Molecule::H2Chain(size_t n, double a, double b, int charge, CartesianDirection axis) {

    return Molecule(NuclearFramework::H2Chain(n, a, b, axis), charge);
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

    os << "Number of electrons: " << molecule.numberOfElectrons() << std::endl;
    os << molecule.nuclearFramework();

    return os;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the sum of all the charges of the nuclei
 */
size_t Molecule::totalNucleicCharge() const {

    return this->nuclear_framework.totalNucleicCharge();
}


/**
 *  @param index1   the index of the first nucleus
 *  @param index2   the index of the second nucleus
 *
 *  @return the distance between the two nuclei at index1 and index2 in bohr
 */
double Molecule::internuclearDistance(const size_t index1, const size_t index2) const {

    return this->nuclearFramework().internuclearDistance(index1, index2);
}


}  // namespace GQCP
