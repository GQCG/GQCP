#ifndef GQCG_MOLECULE_HPP
#define GQCG_MOLECULE_HPP


#include <stdlib.h>
#include <string>
#include <vector>

#include <Atom.hpp>



namespace GQCG {



class Molecule {
private:
    const std::vector<GQCG::Atom> atoms;
    const size_t N;  // number of electrons

    /**
     *  Parse a @param xyz_filename to @return a std::vector<GQCG::Atom>.
     *
     *  The coordinates in the .xyz-file should be in Angstrom: this function converts them immediately to Bohr (a.u.)
     */
    static std::vector<GQCG::Atom> parseXYZFile(const std::string& xyz_filename);


public:
    // CONSTRUCTORS
    /**
     *  Constructor from a given @param xyz_filename
     *      The constructed molecule instance corresponds to a neutral atom (i.e. N = sum of nucleus charges)
     *
     *  IMPORTANT!!! The coordinates of the atoms in the .xyz-file should be in Angstrom, but we convert them internally to Bohr
     */
    explicit Molecule(const std::string& xyz_filename);

    /**
     *  Constructor from a given @param xyz_filename and a @param molecular_charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     *  IMPORTANT!!! The coordinates of the atoms in the .xyz-file should be in Angstrom, but we convert them internally to Bohr
     */
    Molecule(const std::string& xyz_filename, int molecular_charge);

    /**
     *  Constructor from a @param atoms: a given std::vector of GQCG::Atoms
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
     */
    explicit Molecule(const std::vector<GQCG::Atom>& atoms);

    /**
     *  Constructor from a @param atoms: a given std::vector of GQCG::Atoms and a @param molecular_charge
     *      The constructed molecule instance corresponds to an ion:
     *          charge = +1 -> cation (one electron less than the neutral molecule)
     *          charge = 0  -> neutral molecule
     *          charge = -1 -> anion (one electron more than the neutral molecule)
     *
     *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
     */
    Molecule(const std::vector<GQCG::Atom>& atoms, int molecular_charge);


    // OPERATORS
    /**
     *  @return if this is equal to @param other, within the default GQCG::Atom::tolerance_for_comparison for the coordinates of the atoms
     */
    bool operator==(const GQCG::Molecule& other) const;

    /**
     *  Overloading of operator<< for a GQCG::Molecule to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCG::Molecule& molecule);


    // GETTERS
    size_t get_N() const { return this->N; }
    size_t numberOfAtoms() const { return this->atoms.size(); }


    // PUBLIC METHODS
    /**
     *  @return if this is equal to @param other, within the given @param tolerance for the coordinates of the atoms
     */
    bool isEqualTo(const GQCG::Molecule& other, double tolerance=GQCG::Atom::tolerance_for_comparison) const;

    /**
     *  @return the sum of all the charges of the nuclei
     */
    size_t calculateTotalNucleicCharge() const;

    /**
     *  @return the distance between two the two atoms at @param index1 and @param index2
     */
    double calculateInternuclearDistance(size_t index1, size_t index2) const;

    /**
     *  @return the internuclear repulsion energy due to the nuclear framework
     */
    double calculateInternuclearRepulsionEnergy() const;
};



}  // namespace GQCG


#endif  // GQCG_MOLECULE_HPP
