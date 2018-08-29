#ifndef GQCG_MOLECULE_HPP
#define GQCG_MOLECULE_HPP


#include <stdlib.h>
#include <string>
#include <vector>

#include <Atom.hpp>



namespace GQCG {



class Molecule {
private:
    const size_t N;  // number of electrons
    const std::vector<GQCG::Atom> atoms;

    /**
     *  Parse a @param xyz_filename to @return a std::vector<GQCG::Atom>.
     *
     *  The coordinates in the .xyz-file should be in Angstrom: this function converts them immediately to Bohr (a.u.)
     */
    std::vector<GQCG::Atom> parseXYZFile(const std::string& xyz_filename) const;


public:
    // CONSTRUCTORS
    /**
     *  Constructor from a given @param xyz_filename
     *      The constructed molecule instance corresponds to a neutral atom (i.e. N = sum of nucleus charges)
     *
     *  IMPORTANT!!! The coordinates of the atoms in the .xyz-file should be in Angstrom, but we convert them internally to Bohr
     */
    Molecule(std::string xyz_filename);
};



}  // namespace GQCG


#endif  // GQCG_MOLECULE_HPP
