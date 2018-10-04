#include "PlainRHFSCFSolver.hpp"


namespace GQCG {

/*
 *  PRIVATE METHODS
 */

/**
 *  Calculate a new Fock matrix (in AO basis), i.e. this is the 'plain' RHF SCF step
 */
Eigen::MatrixXd PlainRHFSCFSolver::calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) {
    return GQCG::calculateRHFAOFockMatrix(D_AO, this->ham_par);
}


/*
 * CONSTRUCTORS
 */

/**
 *  Constructor based on given @param hamiltonian parameters, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
 */
PlainRHFSCFSolver::PlainRHFSCFSolver(GQCG::HamiltonianParameters ham_par, GQCG::Molecule molecule, double threshold, size_t maximum_number_of_iterations) :
    GQCG::RHFSCFSolver(ham_par, molecule, threshold, maximum_number_of_iterations)
{}


}  // namespace GQCG
