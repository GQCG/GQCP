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
#include "Basis/LibcintInterfacer.hpp"



namespace GQCP {


/*
 *  PUBLIC METHODS - INTERFACING
 */

/**
 *  @param shell_set        the GQCP::ShellSet whose information should be converted
 *
 *  @return the information in a GQCP::ShellSet as a libcint::Container
 */
libcint::RawContainer LibcintInterfacer::convert(const ShellSet& shell_set) const{

    const auto& nuclei = shell_set.nuclei();
    const auto& natm = nuclei.size();
    const auto& nbf = shell_set.numberOfBasisFunctions();
    const auto& nsh = shell_set.size();  // number of shells

    libcint::RawContainer raw_container (natm, nbf, nsh);


    // Configuration of atom-related data
    int offset = libcint::ptr_env_start;  // an offset such that libcint can retrieve the correct index inside the environment, starts at 20

    for (size_t i = 0; i < natm; i++) {
        // Configure a libcint 'atom'
        raw_container.libcint_atm[libcint::charge_of + libcint::atm_slots * i] = static_cast<int>(nuclei[i].charge());  // insert the charge/atomic number
        raw_container.libcint_atm[libcint::ptr_coord + libcint::atm_slots * i] = offset;  // 'pointer' to the coordinates of the atom inside the libcint environment

        // Set the atom-related data into the libcint environment
        raw_container.libcint_env[offset + 0] = nuclei[i].position().x();  // insert the position of the nuclei
        raw_container.libcint_env[offset + 1] = nuclei[i].position().y();
        raw_container.libcint_env[offset + 2] = nuclei[i].position().z();
        offset += 3;
    }



    // Configuration of shell-related data
    int nucleus_index = 0;  // index of the nucleus the shell is centered on
    auto previous_nucleus = shell_set[0].get_nucleus();  // start with the first nucleus
    for (size_t n = 0; n < shell_set.numberOfShells(); n++) {

        auto current_shell = shell_set[n];
        if (current_shell.is_normalized()) {
            throw std::invalid_argument("LibcintInterfacer::convert(const ShellSet&): The libcint integral engine requires a ShellSet with coefficients that do not hold the total normalization factor.");
        }


        // If there's a new nucleus, increment the index
        const auto& current_nucleus = current_shell.get_nucleus();
        if (!Nucleus::equalityComparer()(current_nucleus, previous_nucleus)) {
            nucleus_index++;
            previous_nucleus = current_nucleus;
        }
        raw_container.libcint_bas[libcint::atom_of + libcint::bas_slots * n] = nucleus_index;


        // Set shell-related data into the libcint 'basis' and into the libcint environment
        raw_container.libcint_bas[libcint::ang_of + libcint::bas_slots * n] = static_cast<int>(current_shell.get_l());  // angular momentum
        raw_container.libcint_bas[libcint::nprim_of + libcint::bas_slots * n] = static_cast<int>(current_shell.contractionSize());  // number of primitives
        raw_container.libcint_bas[libcint::nctr_of + libcint::bas_slots * n] = 1;  // apparently, the number of contractions is always 1 (I'm still not sure what the libcint number of contractions means)
        raw_container.libcint_bas[libcint::ptr_exp + libcint::bas_slots * n] = offset;  // pointer to the exponents of the shell inside the libcint environment

        const auto& gaussian_exponents = current_shell.get_gaussian_exponents();
        for (size_t e = 0; e < gaussian_exponents.size(); e++, offset++) {  // also increment offset
            raw_container.libcint_env[offset] = gaussian_exponents[e];
        }


        raw_container.libcint_bas[libcint::ptr_coeff + libcint::bas_slots * n] = offset;  // pointer to the contraction coefficients inside the libcint environment
        // Input NORMALIZED contraction coefficients inside the libcint 'env'
        if (current_shell.are_embedded_normalization_factors_of_primitives()) {
            current_shell.unEmbedNormalizationFactorsOfPrimitives();
        }

        const auto& current_contraction_coefficients = current_shell.get_contraction_coefficients();
        for (size_t c = 0; c < current_contraction_coefficients.size(); c++, offset++) {  // also increment offset
            raw_container.libcint_env[offset] = current_contraction_coefficients[c] * CINTgto_norm(raw_container.libcint_bas[libcint::ang_of + libcint::bas_slots*n], raw_container.libcint_env[raw_container.libcint_bas[libcint::ptr_exp + libcint::bas_slots*n]+c]);  // use libcint to embed the norm of the primitives into the contraction coefficient
        }
    }  // shell loop


    return raw_container;
}


/**
 *  Set the origin for the calculation of all vector-related integrals
 *
 *  @param raw_container        the libcint::RawContainer that holds the data needed by libcint
 *  @param origin               the new origin for the calculation of all vector-related integrals
 */
void LibcintInterfacer::setCommonOrigin(libcint::RawContainer& raw_container, const Vector<double, 3>& origin) const {

    raw_container.libcint_env[libcint::ptr_common_orig + 0] = origin.x();  // input the origin inside the libcint environment
    raw_container.libcint_env[libcint::ptr_common_orig + 1] = origin.y();
    raw_container.libcint_env[libcint::ptr_common_orig + 2] = origin.z();
}



/*
 *  PUBLIC METHODS - INTEGRALS
 */

/**
 *  @param function             the libcint two-electron integral function
 *  @param raw_container        the data libcint needs to perform calculations
 *
 *  @return an array of N OneElectronOperators corresponding to the matrix representations of the N components of the given operator represented by the libcint function
 */
TwoElectronOperator<double> LibcintInterfacer::calculateTwoElectronIntegrals(const Libcint2eFunction& function, libcint::RawContainer& raw_container) const {

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
