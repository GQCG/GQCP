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
#pragma once


#include "Molecule/NuclearFramework.hpp"

#include <cstdlib>


namespace GQCP {


/**
 *  A molecule: a number of electrons and a collection of nuclei arranged in a nuclear framework.
 */
class Molecule {
private:
    size_t N;                            // the number of electrons
    NuclearFramework nuclear_framework;  // the underlying nuclear framework


public:
    // CONSTRUCTORS

    /**
     *  @param nuclear_framework        the nuclear framework that makes up the molecule, with coordinates in bohr
     *  @param charge                   the charge of the molecule:
     *                                      +1 -> cation (one electron less than the neutral molecule)
     *                                       0 -> neutral molecule
     *                                      -1 -> anion (one electron more than the neutral molecule)
     */
    Molecule(const NuclearFramework& nuclear_framework, const int charge = 0);

    /**
     *  @param nuclei       the nuclei that make up the molecule, with coordinates in bohr
     *  @param charge       the charge of the molecule:
     *                          +1 -> cation (one electron less than the neutral molecule)
     *                           0 -> neutral molecule
     *                          -1 -> anion (one electron more than the neutral molecule)
     */
    Molecule(const std::vector<Nucleus>& nuclei, const int charge = 0);


    // NAMED CONSTRUCTORS

    /**
     *  @param n            the number of H nuclei
     *  @param spacing      the internuclear spacing in bohr
     *  @param charge       the total charge
     *
     *  @return a H-chain with equal internuclear spacing
     */
    static Molecule HChain(size_t n, double spacing, int charge = 0, CartesianDirection axis = CartesianDirection::z);

    /**
     *  @param n        the number of H2-molecules
     *  @param a        the internuclear distance in bohr
     *  @param b        the intermolecular distance in bohr
     *  @param charge   the total charge
     *
     *  @return a charged H2-chain
     */
    static Molecule H2Chain(size_t n, double a, double b, int charge = 0, CartesianDirection axis = CartesianDirection::z);

    /**
     *  @param n                the number of hydrogens
     *  @param distance         the distance (in bohr) between neighbouring hydrogen atoms
     * 
     *  @return a regular H-ring where neighbouring hydrogens are separated by the given distance
     */
    static Molecule HRingFromDistance(const size_t n, const double distance, int charge = 0);

    /**
     *  @param n                the number of hydrogens
     *  @param radius           the radius (in bohr) of the circumscribed circle
     * 
     *  @return a regular H-ring whose hydrogens are on the circle with the given radius
     */
    static Molecule HRingFromRadius(const size_t n, const double radius, int charge = 0);

    /**
     *  Construct a molecule based on the content of a given .xyz-file. In an .xyz-file, the molecular coordinates are in Angstrom
     *
     *  @param xyz_filename     the .xyz-file that contains the molecular coordinates in Angstrom
     *  @param charge       the charge of the molecule:
     *                          +1 -> cation (one electron less than the neutral molecule)
     *                           0 -> neutral molecule
     *                          -1 -> anion (one electron more than the neutral molecule)
     */
    static Molecule ReadXYZ(const std::string& xyz_filename, const int charge = 0);


    // OPERATORS

    /**
     *  @param os           the output stream which the molecule should be concatenated to
     *  @param molecule     the molecule that should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);


    // PUBLIC METHODS

    /**
     *  @return the charge of this molecule (in a.u.)
     */
    int charge() const;

    /**
     *  @param index1   the index of the first nucleus
     *  @param index2   the index of the second nucleus
     *
     *  @return the distance between the two nuclei at index1 and index2 in bohr
     */
    double internuclearDistance(const size_t index1, const size_t index2) const { return this->nuclearFramework().internuclearDistance(index1, index2); }

    /**
     *  @return the underlying nuclear framework
     */
    const NuclearFramework& nuclearFramework() const { return this->nuclear_framework; }

    /**
     *  @return the number of atoms in this molecule
     */
    size_t numberOfAtoms() const { return this->nuclear_framework.numberOfNuclei(); }

    /**
     *  @return the number of electron pairs in this molecule. For odd numbers of electrons, the number of electron pairs is equal to that of the (N-1)-even-electron system.
     */
    size_t numberOfElectronPairs() const { return this->numberOfElectrons() / 2; /* this floors to zero */ }

    /**
     *  @return the number of electrons in the molecule
     */
    size_t numberOfElectrons() const { return this->N; }

    /**
     *  @return the sum of all the charges of the nuclei that are in this molecule
     */
    size_t totalNucleicCharge() const { return this->nuclear_framework.totalNucleicCharge(); }
};


}  // namespace GQCP
