#include "PlainRHFSCFSolver.hpp"


namespace GQCG {

/*
 *  PRIVATE METHODS
 */
/**
 *  Calculate a new Fock matrix (expressed in AO basis), i.e. this is the 'plain' RHF SCF step.
 *
 *  The new Fock matrix is calculated as F = H_core + G, i.e. the new Fock matrix is just the AO Fock matrix
 */
Eigen::MatrixXd PlainRHFSCFSolver::calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) {
    return GQCG::calculateRHFAOFockMatrix(D_AO, this->ham_par);
}



/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
 */
PlainRHFSCFSolver::PlainRHFSCFSolver(GQCG::HamiltonianParameters ham_par, GQCG::Molecule molecule, double threshold, size_t maximum_number_of_iterations) :
    GQCG::RHFSCFSolver(ham_par, molecule, threshold, maximum_number_of_iterations)
{}


}  // namespace GQCG
