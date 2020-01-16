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
#include "Processing/Properties/expectation_values.hpp"


namespace GQCP {


/*
 *  MIXED OPERATORS
 */

/**
 *  @param sq_hamiltonian       the electronic Hamiltonian containing one- and two-electron operators
 *  @param one_rdm              the 1-RDM
 *  @param two_rdm              the 2-RDM
 *
 *  @return the expectation value of the given electronic Hamiltonian
 */
double calculateExpectationValue(const SQHamiltonian<double>& sq_hamiltonian, const OneRDM<double>& one_rdm, const TwoRDM<double>& two_rdm) {

    return calculateExpectationValue(sq_hamiltonian.core(), one_rdm)[0] + calculateExpectationValue(sq_hamiltonian.twoElectron(), two_rdm)[0];  // SQHamiltonian contains ScalarSQOperators, so we access with [0]
}


}  // namespace GQCP
