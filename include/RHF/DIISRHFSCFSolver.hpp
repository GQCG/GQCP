//
//  DIISRHFSCFSolver.hpp
//  GQCP
//
//  Created by Laurent Lemmens on 04/10/2018.
//  Copyright © 2018 Ghent Quantum Chemistry Group. All rights reserved.
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
