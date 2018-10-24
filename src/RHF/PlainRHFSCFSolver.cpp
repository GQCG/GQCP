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
#include "RHF/PlainRHFSCFSolver.hpp"


namespace GQCP {

/*
 *  PRIVATE METHODS
 */
/**
 *  Calculate a new Fock matrix (expressed in AO basis), i.e. this is the 'plain' RHF SCF step.
 *
 *  The new Fock matrix is calculated as F = H_core + G, i.e. the new Fock matrix is just the AO Fock matrix
 */
Eigen::MatrixXd PlainRHFSCFSolver::calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) {
    return GQCP::calculateRHFAOFockMatrix(D_AO, this->ham_par);
}



/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
 */
PlainRHFSCFSolver::PlainRHFSCFSolver(GQCP::HamiltonianParameters ham_par, GQCP::Molecule molecule, double threshold, size_t maximum_number_of_iterations) :
    GQCP::RHFSCFSolver(ham_par, molecule, threshold, maximum_number_of_iterations)
{}


}  // namespace GQCP
