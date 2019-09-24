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
 *  @param wavefunction1        a wave function in a product Fock space  
 *  @param wavefunction2        a wave function in a product Fock space containing one fewer electron and the same amount of orbitals that is expressed in the same basis
 *  
 *  @return a vector with the coefficients of a Dyson orbital derived from the difference between the two wave functions expressed in the basis of the wave functions
 */
VectorX<double> calculateDysonAmplitudes(const WaveFunction& wavefunction1, const WaveFunction& wavefunction2) {

    if (wavefunction1.get_fock_space().get_type() != FockSpaceType::ProductFockSpace) {
        throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction1 is not in a product Fock space");
    }

    if (wavefunction2.get_fock_space().get_type() != FockSpaceType::ProductFockSpace) {
        throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction2 is not in a product Fock space");
    }

    const auto fock_space1 = static_cast<const ProductFockSpace&>(wavefunction1.get_fock_space());
    const auto fock_space2 = static_cast<const ProductFockSpace&>(wavefunction2.get_fock_space());

    if (!((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() == 0) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() == 1))) {
        if (!((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() == 1) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() == 0))) {
            throw std::runtime_error("properties::calculateDysonOrbital(WaveFunction, WaveFunction): wavefunction2 is not expressed in a product Fock space with one fewer electron than wavefunction1");
        }
    } 

    const auto fock_space_alpha1 = fock_space1.get_fock_space_alpha();
    const auto fock_space_alpha2 = fock_space2.get_fock_space_alpha();
    const auto fock_space_beta1 = fock_space1.get_fock_space_beta();
    const auto fock_space_beta2 = fock_space2.get_fock_space_beta();

    const auto ci_coeff1 = wavefunction1.get_coefficients();
    const auto ci_coeff2 = wavefunction2.get_coefficients();

    VectorX<double> dyson_coeff = VectorX<double>::Zero(fock_space1.get_K());
    size_t dim_alpha = fock_space_alpha1.get_dimension();
    size_t dim_beta = fock_space_beta1.get_dimension();

    // Given the one electron difference requirement, one of the spin Fock spaces will be identical for both wave functions, identifying that spin Fock space allows for a simple algorithm

    // Beta electron differs by one
    if ((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() == 0) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() == 1)) {

        ONV onv_beta = fock_space_beta1.makeONV(0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings
            int sign = -1;
            for (size_t e_b = 0; e_b < fock_space_beta1.get_N(); e_b++) {  // loop over beta electrons
                sign *= -1;
                size_t p = onv_beta.get_occupation_index(e_b);

                onv_beta.annihilate(p);

                size_t address = fock_space_beta2.getAddress(onv_beta.get_unsigned_representation());

                onv_beta.create(p);

                double coeff = 0;

                for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // alpha Fock space is identical
                    coeff += sign * ci_coeff1(Ia * dim_beta + Ib) * ci_coeff2(address * dim_beta + Ib);
                }

                dyson_coeff(p) += coeff;
            }


            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta1.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop
    }

    // Alpha electrons differs by one
    if ((fock_space1.get_N_alpha() - fock_space2.get_N_alpha() == 1) && (fock_space1.get_N_beta() - fock_space2.get_N_beta() == 0)) {

        ONV onv_alpha = fock_space_alpha1.makeONV(0);

        for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha1.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);

                onv_alpha.annihilate(p);
            
                size_t address = fock_space_alpha2.getAddress(onv_alpha.get_unsigned_representation());

                onv_alpha.create(p);

                double coeff = 0;

                for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // beta Fock space is identical
                    coeff += ci_coeff1(Ia * dim_beta + Ib) * ci_coeff2(address * dim_beta + Ib);
                }

                dyson_coeff(p) += coeff;
            }


            if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
                fock_space_alpha1.setNextONV(onv_alpha);
            }
        }  // alpha address (Ia) loop

    }

    return dyson_coeff;
}


}  // namespace GQCP
