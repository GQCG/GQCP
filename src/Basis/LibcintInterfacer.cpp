#include "Basis/LibcintInterfacer.hpp"



namespace GQCP {



///**
// *  @param libcint_container        the wrapper libcint container
// *
// *  @return the raw information about atoms that libcint can use
// */
//libcint::RawContainer LibcintInterfacer::interface(const libcint::Container& libcint_container) const {
//
//
//
//    const auto& libcint_atoms = libcint_container.atoms;
//    const auto& natm = libcint_atoms.size();
//
//    int offset = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment
//    for (size_t i = 0; i < natm; i++) {
//    }
//
//
//}




/**
 *  @param shell_set        the GQCP::ShellSet whose information should be converted
 *
 *  @return the information in a GQCP::ShellSet as a libcint::Container
 */
libcint::RawContainer LibcintInterfacer::convert(const ShellSet& shell_set) const{

    const auto& atoms = shell_set.atoms();
    const auto& natm = atoms.size();
    const auto& nbf = shell_set.numberOfBasisFunctions();
    const auto& nsh = shell_set.size();  // number of shells

    libcint::RawContainer raw_container (natm, nbf, nsh);


    // Configure 'atm'
    int offset = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment, starts at 20

    for (size_t i = 0; i < natm; i++) {
        raw_container.libcint_atm[CHARGE_OF + ATM_SLOTS * i] = static_cast<int>(atoms[i].atomic_number);  // insert the charge/atomic number
        raw_container.libcint_atm[PTR_COORD + ATM_SLOTS * i] = offset;  // pointer to the coordinates of the atom inside the libcint environment
        raw_container.libcint_env[offset + 0] = atoms[i].position.x();  // insert the position of the atoms
        raw_container.libcint_env[offset + 1] = atoms[i].position.y();
        raw_container.libcint_env[offset + 2] = atoms[i].position.z();
        offset += 3;
    }



    // Configure 'bas'
    int atom_index = 0;  // index of the atom the shell is centered on
    auto previous_atom = shell_set[0].get_atom();

    for (size_t n = 0; n < shell_set.numberOfShells(); n++) {

        auto current_shell = shell_set[n];
        if (current_shell.is_normalized()) {
            throw std::invalid_argument("LibcintInterfacer::convert(const ShellSet&): The libcint integral engine requires a ShellSet with coefficients that do not hold the total normalization factor.");
        }

        const auto& gaussian_exponents = current_shell.get_gaussian_exponents();


        // If there's a new atom, increment the index
        auto current_atom = current_shell.get_atom();
        if (current_atom != previous_atom) {
            atom_index++;
            previous_atom = current_atom;
        }
        raw_container.libcint_bas[ATOM_OF  + BAS_SLOTS * n] = atom_index;


        raw_container.libcint_bas[ANG_OF   + BAS_SLOTS * n] = static_cast<int>(current_shell.get_l());  // angular momentum
        raw_container.libcint_bas[NPRIM_OF + BAS_SLOTS * n] = static_cast<int>(current_shell.contractionSize());  // number of primitives
        raw_container.libcint_bas[NCTR_OF  + BAS_SLOTS * n] = 1;  // number of contractions  // FIXME ????

        raw_container.libcint_bas[PTR_EXP  + BAS_SLOTS * n] = offset;  // pointer to the exponents of the shell inside the libcint environment
        for (size_t e = 0; e < gaussian_exponents.size(); e++, offset++) {  // also increment offset
            raw_container.libcint_env[offset] = gaussian_exponents[e];
        }


        raw_container.libcint_bas[PTR_COEFF + BAS_SLOTS * n] = offset;  // pointer to the contraction coefficients inside the libcint environment


        // Input normalized contraction coefficients inside the libcint 'env'
        if (current_shell.are_embedded_normalization_factors_of_primitives()) {
            current_shell.unEmbedNormalizationFactorsOfPrimitives();
        }

        const auto& current_contraction_coefficients = current_shell.get_contraction_coefficients();


        for (size_t c = 0; c < current_contraction_coefficients.size(); c++, offset++) {  // also increment offset
            raw_container.libcint_env[offset] = current_contraction_coefficients[c] * CINTgto_norm(raw_container.libcint_bas[ANG_OF+BAS_SLOTS*n], raw_container.libcint_env[raw_container.libcint_bas[PTR_EXP+BAS_SLOTS*n]+c]);
        }

    }

    return raw_container;
}




}  // namespace GQCP
