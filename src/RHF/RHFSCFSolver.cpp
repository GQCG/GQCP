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
#include "RHF/RHFSCFSolver.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param ham_par                          the Hamiltonian parameters in AO basis
 *  @param molecule                         the molecule used for the SCF calculation
 *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
 *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
 */
RHFSCFSolver::RHFSCFSolver(HamiltonianParameters ham_par, Molecule molecule, double threshold, size_t maximum_number_of_iterations) :
    ham_par (ham_par),
    molecule (molecule),
    maximum_number_of_iterations (maximum_number_of_iterations),
    threshold (threshold)
{
    // Check if the given molecule has an even number of electrons
    if ((molecule.get_N() % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons.");
    }
}



/*
 *  PUBLIC METHODS
 */
/**
 *  Solve the RHF SCF equations
 */
void RHFSCFSolver::solve() {

    Eigen::MatrixXd H_core = this->ham_par.get_h().get_matrix_representation();
    Eigen::MatrixXd S = this->ham_par.get_S().get_matrix_representation();


    // Obtain an initial guess for the AO density matrix by solving the generalized eigenvalue problem for H_core
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> initial_generalized_eigensolver (H_core, S);
    Eigen::MatrixXd C = initial_generalized_eigensolver.eigenvectors();
    Eigen::MatrixXd D_AO = calculateRHFAO1RDM(C, this->molecule.get_N());


    size_t iteration_counter = 0;
    while (!(this->is_converged)) {
        Eigen::MatrixXd F_AO = this->calculateNewFockMatrix(D_AO);

        // Solve the generalized eigenvalue problem for the Fock matrix to get an improved density matrix
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> generalized_eigensolver (F_AO, S);
        C = generalized_eigensolver.eigenvectors();

        Eigen::MatrixXd D_AO_previous = D_AO;  // store the previous density matrix to be able to check on convergence
        D_AO = calculateRHFAO1RDM(C, this->molecule.get_N());


        // Check for convergence on the AO density matrix
        if ((D_AO - D_AO_previous).norm() <= this->threshold) {
            this->is_converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            OneElectronOperator F (F_AO);
            F.transform(C);  // transform F to the MO basis with C
            if (!(F.get_matrix_representation().isDiagonal())) {
                throw std::runtime_error("The RHF SCF procedure is converged but the MO Fock matrix is not diagonal.");
            }

            // Set the converged solution
            auto electronic_energy = calculateRHFElectronicEnergy(D_AO, H_core, F_AO);
            this->solution = RHF(electronic_energy, C, generalized_eigensolver.eigenvalues());

        } else {  // not converged yet
            iteration_counter++;

            // If we reach more than this->maximum_number_of_iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("The SCF procedure did not converge.");
            }
        }
    }  // while not converged
}


}  // namespace GQCP
