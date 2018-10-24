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
#ifndef DIISRHFSCFSolver_hpp
#define DIISRHFSCFSolver_hpp


#include "RHFSCFSolver.hpp"

#include <deque>



/**
 *  A DIIS RHF SCF solver.
 */
namespace GQCP {

class DIISRHFSCFSolver : public GQCP::RHFSCFSolver {
private:
    const size_t maximum_subspace_dimension;

    std::deque<Eigen::MatrixXd> fock_matrix_deque = {};  // deque of Fock matrices used in the DIIS algorithm
    std::deque<Eigen::MatrixXd> error_matrix_deque = {};  // deque of error matrices used in the DIIS algorithm

    // PRIVATE METHODS
    /**
     *  Calculate a new Fock matrix (in AO basis), i.e. this is the 'DIIS' RHF SCF step.
     */
    Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
     */
    DIISRHFSCFSolver(GQCP::HamiltonianParameters ham_par, GQCP::Molecule molecule, size_t maximum_subspace_dimension = 6, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);
};


}  // namespace GQCP



#endif /* DIISRHFSCFSolver_hpp */
