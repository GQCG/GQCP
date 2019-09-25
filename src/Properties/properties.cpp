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
#include "Properties/properties.hpp"

#include "Properties/expectation_values.hpp"
#include "FockSpace/ProductFockSpace.hpp"


namespace GQCP {


/**
 *  @param dipole_operator      the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm              the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Vector<double, 3> calculateElectronicDipoleMoment(const VectorSQOneElectronOperator<double>& dipole_operator, const OneRDM<double>& one_rdm) {

    auto expectation_values = calculateExpectationValue<3>(dipole_operator, one_rdm);

    Vector<double, 3> electronic_dipole = Eigen::Map<Eigen::Vector3d>(expectation_values.data());
    return electronic_dipole;
}


/**
 *  Calculate the Dyson 'amplitudes' (the coefficients of a Dyson orbital) between two wave function expressed in the same spinor basis 
 * 
 *  @param wavefunction1        a wave function in a product Fock space  
 *  @param wavefunction2        a wave function in a product Fock space containing one fewer electron and the same amount of orbitals that is expressed in the same basis
 *
 *  @return a vector with the Dyson orbital amplitudes  
 */
VectorX<double> calculateDysonAmplitudes(const WaveFunction& wavefunction1, const WaveFunction& wavefunction2) {

    // Check the arguments
    if (wavefunction1.get_fock_space().get_type() != FockSpaceType::ProductFockSpace) {
        throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction1 is not in a product Fock space");
    }

    if (wavefunction2.get_fock_space().get_type() != FockSpaceType::ProductFockSpace) {
        throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction2 is not in a product Fock space");
    }

    const auto& fock_space1 = static_cast<const ProductFockSpace&>(wavefunction1.get_fock_space());
    const auto& fock_space2 = static_cast<const ProductFockSpace&>(wavefunction2.get_fock_space());

    if ((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() != 0) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() != 1)) {
        if ((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() != 1) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() != 0)) {
            throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction2 is not expressed in a product Fock space with one fewer electron than wavefunction1");
        }
    } 

    // Initialize environment variables

    // the passive spin Fock spaces are the Fock spaces that are equal for both wave functions, target Fock spaces have an electron difference of one
    //  we initialize the environments for the case in which they differ one beta electron
    auto passive_fock_space1 = fock_space1.get_fock_space_alpha();
    auto passive_fock_space2 = fock_space2.get_fock_space_alpha();
    auto target_fock_space1 = fock_space1.get_fock_space_beta();
    auto target_fock_space2 = fock_space2.get_fock_space_beta();

    size_t passive_mod1 = target_fock_space1.get_dimension();
    size_t passive_mod2 = target_fock_space2.get_dimension();
    size_t target_mod = 1;

    // If instead the Fock spaces differ by one alpha electron we re-assign the variables to match the algorithm
    if ((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() == 1) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() == 0)) {
        passive_fock_space1 = target_fock_space1;
        passive_fock_space2 = target_fock_space2;
        target_fock_space1 = fock_space1.get_fock_space_alpha();
        target_fock_space2 = fock_space2.get_fock_space_alpha();
        
        passive_mod1 = 1;
        passive_mod2 = 1;
        target_mod = passive_fock_space1.get_dimension();
    }

    const auto& ci_coeffs1 = wavefunction1.get_coefficients();
    const auto& ci_coeffs2 = wavefunction2.get_coefficients();

    VectorX<double> dyson_coeff = VectorX<double>::Zero(fock_space1.get_K());

    // Dyson algorithm
    //  perform a CI iteration over the more occupied Fock space 
    //  for each ONV perform the operation to transform it to a valid ONV for that of the lesser occupied Fock space
    //  based on the address match the correct coefficients
    //  the coefficients products are summed according to the index of the operation on the initial ONV
    ONV onv = target_fock_space1.makeONV(0);

    for (size_t It = 0; It < target_fock_space1.get_dimension(); It++) {  // It loops over addresses of the target Fock space
        int sign = -1;  // operator sign
        for (size_t e = 0; e < target_fock_space1.get_N(); e++) {  // loop over electrons in the onv
            
            sign *= -1; 
            size_t p = onv.get_occupation_index(e);

            // annihilation results in valid onv (N-1) in the target Fock space of wavefunction2
            onv.annihilate(p);

            // retrieve the address of the new onv
            size_t address = target_fock_space2.getAddress(onv.get_unsigned_representation());

            double coeff = 0;
            for (size_t Ip = 0; Ip < passive_fock_space1.get_dimension(); Ip++) {  // passive Fock space is identical and allows for repeat updates
                coeff += sign * ci_coeffs1(It * target_mod + Ip * passive_mod1) * ci_coeffs2(address * target_mod + Ip * passive_mod2);
            }
            dyson_coeff(p) += coeff;

            // restore the onv to allow iteration
            onv.create(p);
        }

        if (It < target_fock_space1.get_dimension() - 1) {  // prevent last permutation to occur
            target_fock_space1.setNextONV(onv);
        }
    }  // target address (It) loop
    
    return dyson_coeff;
}


}  // namespace GQCP
