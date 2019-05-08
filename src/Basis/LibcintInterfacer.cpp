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


/**
 *  @param function             the libcint two-electron integral function
 *  @param raw_container        the data libcint needs to perform calculations
 *
 *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator represented by the libcint function
 */
TwoElectronOperator<double> LibcintInterfacer::calculateTwoElectronIntegrals(const Libint2eFunction& function, libcint::RawContainer& raw_container) const {

    // Initialize the TwoElectronOperator and set to zero
    const auto& nbf = raw_container.nbf;
    TwoElectronOperator<double> g (nbf);
    g.setZero();



    // Calculate the integrals over the shells
    size_t bf1 = 0;  // index of the first basis function in the first shell
    for (size_t sh1 = 0; sh1 < raw_container.nsh; sh1++) {

        int nbf_sh1 = CINTcgto_cart(static_cast<int>(sh1), raw_container.libcint_bas);  // number of basis functions in first shell
        size_t bf2 = 0;  // index of the first basis function in the second shell
        for (size_t sh2 = 0; sh2 < raw_container.nsh; sh2++) {

            int nbf_sh2 = CINTcgto_cart(static_cast<int>(sh2), raw_container.libcint_bas);  // number of basis functions in first shell
            size_t bf3 = 0;  // index of the first basis function in the third shell
            for (size_t sh3 = 0; sh3 < raw_container.nsh; sh3++) {

                int nbf_sh3 = CINTcgto_cart(static_cast<int>(sh3), raw_container.libcint_bas);  // number of basis functions in first shell
                size_t bf4 = 0;  // index of the first basis function in the fourth shell
                for (size_t sh4 = 0; sh4 < raw_container.nsh; sh4++) {

                    int nbf_sh4 = CINTcgto_cart(static_cast<int>(sh4), raw_container.libcint_bas);  // number of basis functions in fourth shell



                    int shell_indices[4];  // indices of the shells to be calculated over
                    shell_indices[0] = static_cast<int>(sh1);
                    shell_indices[1] = static_cast<int>(sh2);
                    shell_indices[2] = static_cast<int>(sh3);
                    shell_indices[3] = static_cast<int>(sh4);


                    double buf[nbf_sh1 * nbf_sh2 * nbf_sh3 * nbf_sh4];  // buffer where the integrals are calculated to
                    function(buf, shell_indices, raw_container.libcint_atm, raw_container.natm, raw_container.libcint_bas, raw_container.nbf, raw_container.libcint_env, nullptr);  // TODO: what is result is zero: skip a shell?


                    for (size_t f1 = 0; f1 < nbf_sh1; f1++) {
                        for (size_t f2 = 0; f2 < nbf_sh2; f2++) {
                            for (size_t f3 = 0; f3 < nbf_sh3; f3++) {
                                for (size_t f4 = 0; f4 < nbf_sh4; f4++) {
                                    const double& computed_integral = buf[f1 + nbf_sh1 * (f2 + nbf_sh2 * (f3 + nbf_sh3 * f4))];  // integrals are packed in column-major form

                                    // Two-electron integrals are given in CHEMIST'S notation: (11|22)
                                    g(f1 + bf1, f2 + bf2, f3 + bf3, f4 + bf4) = computed_integral;
                                }
                            }
                        }
                    } // data access loops

                    bf4 += nbf_sh4;  // update the 'first basis function' with the encountered number of basis functions
                }

                bf3 += nbf_sh3;  // update the 'first basis function' with the encountered number of basis functions
            }

            bf2 += nbf_sh2;  // update the 'first basis function' with the encountered number of basis functions
        }

        bf1 += nbf_sh1;  // update the 'first basis function' with the encountered number of basis functions
    }  // shell loops


    return g;
}


}  // namespace GQCP
