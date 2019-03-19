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
 *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
 *  @param molecule             the molecule containing the atoms on which the shells should be centered
 */
ShellSet::ShellSet(const std::string& basisset_name, const Molecule& molecule) {

    // Since we haven't implemented our own BasisSet class, we use libint to read in the basisset specification file and create the shells
    // TODO no longer use libint2 to read this

    const auto& atoms = molecule.get_atoms();

    libint2::BasisSet libint_basis (basisset_name, LibintCommunicator::get().interface(atoms));

    // copy Libint2 shells to GQCP::Shells
    this->reserve(libint_basis.size());
    for (const auto& libint_shell : libint_basis) {




        std::vector<Contraction> contractions;
        contractions.reserve(libint_shell.contr.size());

        // copy libint2 contractions (contr) to gqcp contractions
        for (const auto &contraction : libint_shell.contr) {
            contractions.push_back({static_cast<size_t>(contraction.l), contraction.coeff});
        }

        // Libint2 only stores the origin of the shell, so we have to find the atom corresponding to the copied shell's origin
        Atom corresponding_atom;
        for (const Atom &atom : atoms) {
            Eigen::Map<const Eigen::Matrix<double, 3, 1>> libint_origin_map(libint_shell.O.data());
            if (atom.position.isApprox(libint_origin_map)) {
                corresponding_atom = atom;
                break;
            }
        }

        this->emplace_back(corresponding_atom, libint_shell.alpha, contractions);
    }
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of shells in this basisset
 */
size_t ShellSet::numberOfShells() const {
    return this->size();
}


/**
 *  @return the number of basis functions in this basisset
 */
size_t ShellSet::numberOfBasisFunctions() const {
    size_t number_of_basis_functions = 0;
    for (const auto& shell : *this) {
        for (const auto& contraction : shell.get_contractions()) {
            number_of_basis_functions += contraction.numberOfBasisFunctions();
        }
    }
    return number_of_basis_functions;
}


/**
 *  @return an ordered vector of the unique atoms in this basisset
 */
std::vector<Atom> ShellSet::atoms() const {

    std::vector<Atom> atoms {};

    // Append every unique atom in this basisset's shells
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
