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
#ifndef PlainRHFSCFSolver_hpp
#define PlainRHFSCFSolver_hpp


#include "RHFSCFSolver.hpp"


namespace GQCP {


/**
 *  A plain RHF SCF solver.
 */
class PlainRHFSCFSolver : public GQCP::RHFSCFSolver {
private:
    /**
     *  Calculate a new Fock matrix (expressed in AO basis), i.e. this is the 'plain' RHF SCF step.
     *
     *  The new Fock matrix is calculated as F = H_core + G, i.e. the new Fock matrix is just the AO Fock matrix
     */
    Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
     */
    PlainRHFSCFSolver(GQCP::HamiltonianParameters ham_par, GQCP::Molecule molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);
};


}  // namespace GQCP


#endif /* PlainRHFSCFSolver_hpp */
