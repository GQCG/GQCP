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
#ifndef expectation_values_hpp
#define expectation_values_hpp

#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"


namespace GQCP {


/*
 *  ONE-ELECTRON OPERATORS
 */

/**
 *  @param one_op       the one-electron operator whose expectation value should be calculated
 *  @param one_rdm      the 1-RDM that represents the wave function
 *
 *  @return the expectation value of the one-electron operator, with the given 1-RDM
 */
double calculateExpectationValue(const GQCP::OneElectronOperator& one_op, const GQCP::OneRDM& one_rdm);

/**
 *  @tparam N           the number of components of the one-electron operator
 *
 *  @param one_ops      the components of the one-electron operator
 *  @param one_rdm      the 1-RDM that represents the wave function
 *
 *  @return the expectation values of all components of the one-electron operator
 */
template <size_t N>
std::array<double, N> calculateExpectationValues(const std::array<GQCP::OneElectronOperator, N>& one_ops, const GQCP::OneRDM& one_rdm);



/*
 *  TWO-ELECTRON OPERATORS
 */

/**
 *  @param two_op       the two-electron operator whose expectation value should be calculated
 *  @param two_rdm      the 2-RDM that represents the wave function
 *
 *  @return the expectation value of the one-electron operator, with the given 2-RDM: this includes the prefactor 1/2
 */
double calculateExpectationValue(const GQCP::TwoElectronOperator& two_op, const GQCP::TwoRDM& two_rdm);



/*
 *  MIXED OPERATORS
 */

/**
 *  @param ham_par      the Hamiltonian parameters containing the one- and two-electron integrals
 *  @param one_rdm      the 1-RDM
 *  @param two_rdm      the 2-RDM
 *
 *  @return the expectation value of the 'Hamiltonian' represented by the Hamiltonian parameters
 */
double calculateExpectationValue(const GQCP::HamiltonianParameters& ham_par, const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm);


}  // namespace GQCP



#endif /* expectation_values_hpp */
