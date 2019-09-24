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
#pragma once


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"


namespace GQCP {


/*
 *  ONE-ELECTRON OPERATORS
 */

/**
 *  @tparam Components      the number of components of the one-electron operator
 *
 *  @param one_op           the one-electron operator
 *  @param one_rdm          the 1-RDM that represents the wave function
 *
 *  @return the expectation values of all components of the one-electron operator
 */
template <size_t Components>
std::array<double, Components> calculateExpectationValue(const SQOneElectronOperator<double, Components>& one_op, const OneRDM<double>& one_rdm) {


    if (one_op.dimension() != one_rdm.cols()) {
        throw std::invalid_argument("calculateExpectationValue(ScalarSQOneElectronOperator<double>, OneRDM<double>): The given one-electron integrals are not compatible with the 1-RDM.");
    }


    std::array<double, Components> expectation_values {};  // zero initialization
    for (size_t i = 0; i < Components; i++) {
        expectation_values[i] = (one_op.parameters(i) * one_rdm).trace();
    }

    return expectation_values;
}


/**
 *  @tparam Components      the number of components of the two-electron operator
 * 
 *  @param two_op           the two-electron operator
 *  @param two_rdm          the 2-RDM that represents the wave function
 *
 *  @return the expectation values of all the components of the two-electron operator, with the given 2-RDM: this includes the prefactor 1/2
 */
template <size_t Components>
std::array<double, Components> calculateExpectationValue(const SQTwoElectronOperator<double, Components>& two_op, const TwoRDM<double>& two_rdm) {

    if (two_op.dimension() != two_rdm.dimension()) {
        throw std::invalid_argument("calculateExpectationValue(ScalarSQTwoElectronOperator<double>, TwoRDM<double>): The given two-electron integrals are not compatible with the 2-RDM.");
    }


    std::array<double, Components> expectation_values {};  // zero initialization
    for (size_t i = 0; i < Components; i++) {

        // Specify the contractions for the relevant contraction of the two-electron integrals and the 2-RDM
        //      0.5 g(p q r s) d(p q r s)
        Eigen::array<Eigen::IndexPair<int>, 4> contractions = {Eigen::IndexPair<int>(0,0), Eigen::IndexPair<int>(1,1), Eigen::IndexPair<int>(2,2), Eigen::IndexPair<int>(3,3)};
        //      Perform the contraction
        Eigen::Tensor<double, 0> contraction = 0.5 * two_op.parameters(i).contract(two_rdm.Eigen(), contractions);

        // As the contraction is a scalar (a tensor of rank 0), we should access by (0).
        expectation_values[i] = contraction(0);
    }

    return expectation_values;
}



/*
 *  MIXED OPERATORS
 */

/**
 *  @param sq_hamiltonian       the Hamiltonian in an orthonormal basis
 *  @param one_rdm              the 1-RDM
 *  @param two_rdm              the 2-RDM
 *
 *  @return the expectation value of the Hamiltonian
 */
double calculateExpectationValue(const SQHamiltonian<double>& sq_hamiltonian, const OneRDM<double>& one_rdm, const TwoRDM<double>& two_rdm);


}  // namespace GQCP
