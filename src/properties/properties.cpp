// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#include "properties/properties.hpp"

#include "properties/expectation_values.hpp"


namespace GQCP {


/**
 *  @param dipole_operator      the three components of the Cartesian dipole integrals in the orthonormal basis in which the 1-RDM is expressed
 *  @param one_rdm              the 1-RDM
 *
 *  @return the three Cartesian components of the electronic electric dipole moment
 */
Eigen::Vector3d calculateElectronicDipoleMoment(const std::array<GQCP::OneElectronOperator, 3>& dipole_operator, const GQCP::OneRDM& one_rdm) {

    auto expectation_values = calculateExpectationValues<3>(dipole_operator, one_rdm);

    Eigen::Vector3d electronic_dipole = Eigen::Map<Eigen::Vector3d>(expectation_values.data());
    return electronic_dipole;
}


/**
 *  Calculates the Mulliken operator for Hamiltonian parameters and a set of GTOs indexes
 *
 *  @param ham_par      the Hamiltonian parameters
 *  @param gto_list     indexes of the original GTOs on which the Mulliken populations are dependant
 *
 *  @return the Mulliken operator for a set of GTOs
 */
OneElectronOperator calculateMullikenOperator(const GQCP::HamiltonianParameters& ham_par, const Vectoru& gto_list) {

    if (!ham_par.get_ao_basis()) {
        throw std::invalid_argument("The Hamiltonian parameters has no underlying GTO basis, Mulliken analysis is not possible.");
    }

    if (gto_list.size() > ham_par.get_K()) {
        throw std::invalid_argument("To many GTOs are selected");
    }

    Eigen::MatrixXd p_a = Eigen::MatrixXd::Zero(ham_par.get_K(), ham_par.get_K());

    for (size_t index : gto_list) {
        if (gto_list.size() >= ham_par.get_K()) {
            throw std::invalid_argument("GTO index is too large");
        }

        p_a(index, index) = 1;
    }

    Eigen::MatrixXd C_inverse = ham_par.get_C().inverse();

    //  Formula for the Mulliken matrix
    Eigen::MatrixXd mulliken_matrix = (C_inverse * p_a * ham_par.get_C() + ham_par.get_C() * p_a * C_inverse) / 2;

    return OneElectronOperator(mulliken_matrix);

}

}  // namespace GQCP
