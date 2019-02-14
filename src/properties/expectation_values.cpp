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
#include "properties/expectation_values.hpp"


namespace GQCP {


/*
 *  ONE-ELECTRON OPERATORS
 */

/**
 *  @param one_op       the one-electron operator whose expectation value should be calculated
 *  @param one_rdm      the 1-RDM that represents the wave function
 *
 *  @return the expectation value of the one-electron operator
 */
double calculateExpectationValue(const OneElectronOperator<double>& one_op, const OneRDM& one_rdm) {

    if (one_op.get_dim() != one_rdm.get_dim()) {
        throw std::invalid_argument("The given one-electron integrals are not compatible with the 1-RDM.");
    }

//    auto h = one_op.get_matrix_representation();
    auto D = one_rdm.get_matrix_representation();

    return (one_op * D).trace();
}



/*
 *  TWO-ELECTRON OPERATORS
 */

/**
 *  @param two_op       the two-electron operator whose expectation value should be calculated
 *  @param two_rdm      the 2-RDM that represents the wave function
 *
 *  @return the expectation value of the one-electron operator: this includes the prefactor 1/2
 */
double calculateExpectationValue(const TwoElectronOperator<double>& two_op, const TwoRDM& two_rdm) {

    if (two_op.get_dim() != two_rdm.get_dim()) {
        throw std::invalid_argument("The given two-electron integrals are not compatible with the 2-RDM.");
    }

//    auto g = two_op.get_matrix_representation();
    auto d = two_rdm.get_matrix_representation();

    // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
    //      0.5 g(p q r s) d(p q r s)
    Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
    //      Perform the contraction
    Eigen::Tensor<double, 0> contraction = 0.5 * two_op.contract(d, contractions);

    // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
    return contraction(0);
}



/*
 *  MIXED OPERATORS
 */

/**
 *  @param ham_par      the Hamiltonian parameters containing the scalar interaction term and the one- and two-electron integrals
 *  @param one_rdm      the 1-RDM
 *  @param two_rdm      the 2-RDM
 *
 *  @return the expectation value of the 'Hamiltonian' represented by the Hamiltonian parameters
 */
double calculateExpectationValue(const HamiltonianParameters<double>& ham_par, const OneRDM& one_rdm, const TwoRDM& two_rdm) {

    return ham_par.get_scalar() + calculateExpectationValue(ham_par.get_h(), one_rdm) + calculateExpectationValue(ham_par.get_g(), two_rdm);
}


}  // namespace GQCP
