#include "Molecule.hpp"

#include <fstream>
#include <sstream>

#include "elements.hpp"
#include "units.hpp"


namespace GQCG {


/*
 *  PRIVATE METHODS
 */

/**
 *  Parse a @param xyz_filename to @return a std::vector<GQCG::Atom>.
 *
 *  The coordinates in the .xyz-file should be in Angstrom: this function converts them immediately to Bohr (a.u.)
 */
std::vector<GQCG::Atom> Molecule::parseXYZFile(const std::string& xyz_filename) {

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
        std::vector<GQCG::Atom> atoms;
        while (std::getline(input_file_stream, line)) {
            std::string symbol;
            double x_angstrom, y_angstrom, z_angstrom;

            std::istringstream iss (line);
            iss >> symbol >> x_angstrom >> y_angstrom >> z_angstrom;

            // Convert the (x,y,z)-coordinates that are in Angstrom to Bohr
            double x_bohr = GQCG::units::angstrom_to_bohr(x_angstrom);
            double y_bohr = GQCG::units::angstrom_to_bohr(y_angstrom);
            double z_bohr = GQCG::units::angstrom_to_bohr(z_angstrom);

            atoms.emplace_back(GQCG::elements::element_to_atomic_number(symbol), x_bohr, y_bohr, z_bohr);
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
 *  Constructor from a @param atoms: a given std::vector of GQCG::Atoms and a @param molecular_charge
 *      The constructed molecule instance corresponds to an ion:
 *          charge = +1 -> cation (one electron less than the neutral molecule)
 *          charge = 0  -> neutral molecule
 *          charge = -1 -> anion (one electron more than the neutral molecule)
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<GQCG::Atom>& atoms, int molecular_charge) :
    atoms (atoms),
    N (this->calculateTotalNucleicCharge() - molecular_charge)
{
    if (molecular_charge > 0) {
        if (this->calculateTotalNucleicCharge() < molecular_charge) {
            throw std::invalid_argument("You cannot create a molecule with these atoms and this much of a total positive charge.");
        }
    }
}


/**
 *  Constructor from a @param atoms: a given std::vector of GQCG::Atoms
 *
 *  IMPORTANT!!! The coordinates of the atoms should be input in Bohr.
 */
Molecule::Molecule(const std::vector<GQCG::Atom>& atoms) :
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
 *  PUBLIC METHODS
 */
/**
 * @return the sum of all the charges of the nuclei
 */
size_t Molecule::calculateTotalNucleicCharge() const {

    size_t nucleic_charge = 0;

    for (const auto& atom : this->atoms) {
        nucleic_charge += atom.atomic_number;
    }

    return nucleic_charge;
}


}  // namespace GQCG