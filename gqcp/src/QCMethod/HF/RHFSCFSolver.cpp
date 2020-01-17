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
#include "QCMethod/HF/BaseRHFSCFSolver.hpp"

#include "QCModel/RHF/RHF.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param sq_hamiltonian                   the Hamiltonian parameters in AO basis
 *  @param spinor_basis                     the spinor basis
 *  @param molecule                         the molecule used for the SCF calculation
 *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
 *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
 */
BaseRHFSCFSolver::BaseRHFSCFSolver(const SQHamiltonian<double>& sq_hamiltonian, const RSpinorBasis<double, GTOShell>& spinor_basis, const Molecule& molecule, double threshold, size_t maximum_number_of_iterations) :
    sq_hamiltonian (sq_hamiltonian),
    spinor_basis (spinor_basis),
    molecule (molecule),
    maximum_number_of_iterations (maximum_number_of_iterations),
    threshold (threshold)
{
    // Check if the given molecule has an even number of electrons
    if ((molecule.numberOfElectrons() % 2) != 0) {
        throw std::invalid_argument("BaseRHFSCFSolver::BaseRHFSCFSolver(): The given molecule has an odd number of electrons.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  Solve the RHF SCF equations
 */
void BaseRHFSCFSolver::solve() {

    const auto& H_core = this->sq_hamiltonian.core().parameters();
    const auto S = this->spinor_basis.overlap().parameters();


    // Obtain an initial guess for the AO density matrix by solving the generalized eigenvalue problem for H_core
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> initial_generalized_eigensolver (H_core, S);
    TransformationMatrix<double> C_initial = initial_generalized_eigensolver.eigenvectors();

    this->solve(C_initial);
}


/**
 *  Solve the RHF SCF equations using an initial guess
 * 
 *  @param C        the initial guess for the canonical RHF coefficient matrix
 */
void BaseRHFSCFSolver::solve(const TransformationMatrix<double>& C_initial) {

    const auto& H_core = this->sq_hamiltonian.core();
    const auto S = this->spinor_basis.overlap().parameters();

    auto C = C_initial;
    auto D_AO = QCModel::RHF<double>::calculateScalarBasis1RDM(C, this->molecule.numberOfElectrons());


    size_t iteration_counter = 0;
    while (!(this->is_converged)) {
        auto F_AO = this->calculateNewFockMatrix(D_AO);

        // Solve the generalized eigenvalue problem for the Fock matrix to get an improved density matrix
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> generalized_eigensolver (F_AO.parameters(), S);
        C = generalized_eigensolver.eigenvectors();

        OneRDM<double> D_AO_previous = D_AO;  // store the previous density matrix to be able to check on convergence
        D_AO = QCModel::RHF<double>::calculateScalarBasis1RDM(C, this->molecule.numberOfElectrons());


        // Check for convergence on the AO density matrix
        if ((D_AO - D_AO_previous).norm() <= this->threshold) {
            this->is_converged = true;

            // After the SCF procedure, we end up with canonical spatial orbitals, i.e. the Fock matrix should be diagonal in this basis
            ScalarSQOneElectronOperator<double> F = F_AO;
            F.transform(C);  // transform F to the MO basis with C
            if (!(F.parameters().isDiagonal())) {
                throw std::runtime_error("BaseRHFSCFSolver::solve(): The RHF SCF procedure is converged but the MO Fock matrix is not diagonal.");
            }

            // Set the converged solution
            auto electronic_energy = QCModel::RHF<double>::calculateElectronicEnergy(D_AO, H_core, F_AO);
            this->solution = RHF(electronic_energy, TransformationMatrix<double>(C), generalized_eigensolver.eigenvalues());

        } else {  // not converged yet
            iteration_counter++;

            // If we reach more than this->maximum_number_of_iterations, the system is considered not to be converging
            if (iteration_counter >= this->maximum_number_of_iterations) {
                throw std::runtime_error("BaseRHFSCFSolver::solve(): The SCF procedure did not converge.");
            }
        }
    }  // while not converged
}


}  // namespace GQCP
