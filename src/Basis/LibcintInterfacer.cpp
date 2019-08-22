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
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"


namespace GQCP {


/*
 *  PUBLIC METHODS - INTERFACING
 */

/**
 *  @param shell_set        the GQCP::ShellSet whose information should be converted
 *
 *  @return the information in a GQCP::ShellSet as a libcint::Container
 */
libcint::RawContainer LibcintInterfacer::convert(const ShellSet<GTOShell>& shell_set) const {

    const auto& nuclei = shell_set.nuclei();
    const auto& natm = nuclei.size();
    const auto& nbf = shell_set.numberOfBasisFunctions();
    const auto& nsh = shell_set.numberOfShells();

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
    const auto& shellset_vector = shell_set.asVector();
    auto previous_nucleus = shellset_vector[0].get_nucleus();  // start with the first nucleus
    for (size_t n = 0; n < shell_set.numberOfShells(); n++) {

        auto current_shell = shellset_vector[n];
        if (current_shell.is_normalized()) {
            throw std::invalid_argument("LibcintInterfacer::convert(const ShellSet<GTOShell>&): The libcint integral engine requires a ShellSet with coefficients that do not hold the total normalization factor.");
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
 *  PUBLIC METHODS - INTEGRAL FUNCTIONS
 */

/**
 *  @param op           the overlap operator
 * 
 *  @return the Libcint one-electron function that corresponds to the overlap operator
 */
Libcint1eFunction LibcintInterfacer::oneElectronFunction(const OverlapOperator& op) const {
    return cint1e_ovlp_cart;
}


/**
 *  @param op               the kinetic operator
 * 
 *  @return the Libcint one-electron function that corresponds to the kinetic operator
 */
Libcint1eFunction LibcintInterfacer::oneElectronFunction(const KineticOperator& op) const {
    return cint1e_kin_cart;
}


/**
 *  @param op               the nuclear attraction operator
 * 
 *  @return the Libcint one-electron function that corresponds to the nuclear attraction operator
 */
Libcint1eFunction LibcintInterfacer::oneElectronFunction(const NuclearAttractionOperator& op) const {
    return cint1e_nuc_cart;
}


/**
 *  @param op               the electronic electric dipole operator
 * 
 *  @return the Libcint one-electron function that corresponds to the electronic electric dipole operator
 */
Libcint1eFunction LibcintInterfacer::oneElectronFunction(const ElectronicDipoleOperator& op) const {
    return cint1e_r_cart;
}


/**
 *  @param op               the Coulomb repulsion operator
 * 
 *  @return the Libcint two-electron function that corresponds to the Coulomb repulsion dipole operator
 */
Libcint2eFunction LibcintInterfacer::twoElectronFunction(const CoulombRepulsionOperator& op) const {

    return cint2e_cart;
}


/**
 *  @param op               the Coulomb repulsion operator
 * 
 *  @return the Libcint two-electron optimizer function that corresponds to the Coulomb repulsion dipole operator
 */
Libcint2eOptimizerFunction LibcintInterfacer::twoElectronOptimizerFunction(const CoulombRepulsionOperator& op) const {

    return cint2e_cart_optimizer;
}


}  // namespace GQCP
