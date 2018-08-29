#include "Molecule.hpp"

#include <fstream>
#include <sstream>

#include "elements.hpp"
#include "units.hpp"


namespace GQCG {


/**
 *  Parse a @param xyz_filename to @return a std::vector<GQCG::Atom>.
 *
 *  The coordinates in the .xyz-file should be in Angstrom: this function converts them immediately to Bohr (a.u.)
 */
std::vector<GQCG::Atom> Molecule::parseXYZFile(const std::string& xyz_filename) const {

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


}  // namespace GQCG