//
//  DIISRHFSCFSolver.cpp
//  gqcg
//
//  Created by Laurent Lemmens on 04/10/2018.
//  Copyright © 2018 Ghent Quantum Chemistry Group. All rights reserved.
//

#include "RHF/DIISRHFSCFSolver.hpp"


namespace GQCG {


/*
 *  PRIVATE METHODS
 */

/**
 *  Calculate a new Fock matrix (in AO basis), i.e. this is the 'DIIS' RHF SCF step.
 */
Eigen::MatrixXd DIISRHFSCFSolver::calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) {

    Eigen::MatrixXd S = this->ham_par.S.get_matrix_representation();

    // Calculate the Fock matrix based off the density matrix
    auto f_AO = GQCG::calculateRHFAOFockMatrix(D_AO, this->ham_par);


    // Update deques for the DIIS procedure
    this->fock_matrix_deque.emplace_back(f_AO);
    Eigen::MatrixXd error_matrix = f_AO * D_AO * S - S * D_AO * f_AO;
    this->error_matrix_deque.emplace_back(error_matrix);


    // Do DIIS when the current subspace dimension is large enough
    // Collapse the subspace, if needed
    size_t n = error_matrix_deque.size();  // n is the current subspace dimension
    if (n == this->maximum_subspace_dimension) {

        // Initialize and calculate the augmented B matrix
        Eigen::MatrixXd B = -1 * Eigen::MatrixXd::Ones(n+1,n+1);  // +1 for the multiplier
        B(n,n) = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                // B(i,j) = Tr(e_i^T e_j)
                B(i,j) = (this->error_matrix_deque[i].transpose() * this->error_matrix_deque[j]).trace();
            }
        }

        // Initialize the RHS of the system of equations
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n+1);  // +1 for the multiplier
        b(n) = -1;  // the last entry of b is accessed through n: dimension of b is n+1 - 1 because of computers


        // Solve the DIIS non-linear equations
        Eigen::VectorXd y = B.inverse() * b;


        // Use the coefficients that are in y to construct 'the best' Fock matrix as a linear combination of previous Fock matrices
        f_AO = Eigen::MatrixXd::Zero(S.cols(), S.cols());
        for (size_t i = 0; i < n; i++) {  // n is the current subspace dimension (not equal to the size of the augmented B matrix)
            f_AO += y(i) * this->fock_matrix_deque[i];
        }

        // Remove the oldest entries, which means that we collapse every iteration once the dimension is large enough
        this->fock_matrix_deque.pop_front();
        this->error_matrix_deque.pop_front();
    }  // subspace collapse

    return f_AO;
}



/*
 *  CONSTRUCTORS
 */
/**
 *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
 */
DIISRHFSCFSolver::DIISRHFSCFSolver(GQCG::HamiltonianParameters ham_par, GQCG::Molecule molecule, size_t maximum_subspace_dimension, double threshold, size_t maximum_number_of_iterations) :
    RHFSCFSolver(ham_par, molecule, threshold, maximum_number_of_iterations),
    maximum_subspace_dimension (maximum_subspace_dimension)
{}


}  // namespace GQCG
