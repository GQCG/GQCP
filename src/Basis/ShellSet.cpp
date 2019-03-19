#include "Basis/ShellSet.hpp"

#include "LibintCommunicator.hpp"

#include <algorithm>


namespace GQCP {



/*
 *  CONSTRUCTORS
 */

/**
 *  Construct a ShellSet by placing the shells corresponding to the basisset information on every atom of the molecule
 *
 *  @param molecule             the molecule containing the atoms on which the shells should be centered
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 */
ShellSet::ShellSet(const Molecule& molecule, const std::string& basisset_name) {

    // Since we haven't implemented our own BasisSet class, we use libint to read in the basisset specification file and create the shells
    // TODO no longer use libint2 to read this

    const auto& atoms = molecule.get_atoms();

    libint2::BasisSet libint_basis (basisset_name, LibintCommunicator::get().interface(atoms));

    *this = LibintCommunicator::get().interface(libint_basis, atoms);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of shells in this shell set
 */
size_t ShellSet::numberOfShells() const {
    return this->size();
}


/**
 *  @return the number of basis functions in this shell set
 */
size_t ShellSet::numberOfBasisFunctions() const {

    size_t nbf {};
    for (const auto& shell : *this) {
        nbf += shell.numberOfBasisFunctions();
    }
    return nbf;
}


/**
 *  @return an ordered vector of the unique atoms in this shell set
 */
std::vector<Atom> ShellSet::atoms() const {

    std::vector<Atom> atoms {};

    // Append every unique atom in this shell set's shells
    for (const auto& shell : *this) {
        auto atom = shell.get_atom();

        auto p = std::find(atoms.begin(), atoms.end(), atom);
        if (p == atoms.end()) {  // if unique
            atoms.push_back(atom);
        }
    }

    return atoms;
}


/**
 *  @param shell_index      the index of the shell
 *
 *  @return the (total basis function) index that corresponds to the first basis function in the given shell
 */
size_t ShellSet::basisFunctionIndex(size_t shell_index) const {

    size_t bf_index {};

    // Count the number of basis functions before
    for (size_t i = 0; i < shell_index; i++) {
        bf_index += this->operator[](i).numberOfBasisFunctions();
    }

    return bf_index;
}


}  // namespace GQCP
