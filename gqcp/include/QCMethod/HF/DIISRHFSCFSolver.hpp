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
#pragma once


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Optimization/IterativeSolver.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "QCModel/HF/RHF.hpp"
#include "QCMethod/HF/BaseRHFSCFSolver.hpp"

#include <Eigen/Dense>


namespace GQCP {

/**
 *  A restricted Hartree-Fock self-consistent field (RHF SCF) solver that uses a DIIS accelerator on the Fock matrix.
 * 
 *  This is not a regular DIIS solver, as it does not construct a linear combination of the iterates (i.e. the coefficient matrices): it uses a linear combination of previous Fock matrices, from which the next iterate is calculated.
 * 
 *  @tparam _ExpansionScalar        the type of scalar that is used to describe the expansion coefficients
 */
template <typename _ExpansionScalar>
class DIISRHFSCFSolver :
    public DIISSolver<TransformationMatrix<_ExpansionScalar>, SquareMatrix<_ExpansionScalar>, _ExpansionScalar>,
    public BaseRHFSCFSolver<_ExpansionScalar> {

public:
    using ExpansionScalar = _ExpansionScalar;
    using BaseDIISSolver = DIISSolver<TransformationMatrix<ExpansionScalar>, SquareMatrix<ExpansionScalar>, ExpansionScalar>;
    using BaseRHFSCFSolver = BaseRHFSCFSolver<ExpansionScalar>;


private:
    double threshold;  // the convergence threshold on the norm of the difference of consecutive density matrices

    OneRDM<ExpansionScalar> D_previous;  // expressed in the scalar orbital basis
    OneRDM<ExpansionScalar> D_current;  // expressed in the scalar orbital basis

    std::deque<ScalarSQOneElectronOperator<ExpansionScalar>> fock_deque;  // expressed in the scalar orbital basis


public:

    /*
     * CONSTRUCTORS
     */

    /**
     *  A general constructor from the properties
     * 
     *  @param C_initial                            the initial guess for the coefficient matrix
     *  @param N                                    the total number of electrons
     *  @param S                                    the overlap matrix of the spinor basis
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform   
     *  @param minimum_subspace_dimension           the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
     *  @param threshold                            the convergence threshold on the norm of the difference of consecutive density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     */
    DIISRHFSCFSolver(const TransformationMatrix<ExpansionScalar>& C_initial, const size_t N, const SquareMatrix<ExpansionScalar>& S, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian, const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6, const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) :
        BaseDIISSolver(C_initial, maximum_number_of_iterations, minimum_subspace_dimension, maximum_subspace_dimension),
        BaseRHFSCFSolver(N, S, sq_hamiltonian),
        threshold (threshold)
    {
        // Given the initial coefficient matrix, we can calculate the initial density matrix
        this->D_previous = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(this->iterate, this->N);
        this->D_current = this->D_previous;  // at the start of the algorithm, they are equal
    }


    /*
     *  NAMED CONSTRUCTORS
     */

    /**
     *  @param spinor_basis                         the spinor basis, containing the (non-orthogonal) underlying scalar orbitals
     *  @param sq_hamiltonian                       the Hamiltonian expressed in that spinor basis
     *  @param N                                    the total number of electrons
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform   
     *  @param minimum_subspace_dimension           the minimum number of iterates that have to be in the subspace before enabling the DIIS acceleration
     *  @param threshold                            the convergence threshold on the norm of the difference of consecutive density matrices
     *  @param maximum_number_of_iterations         the maximum number of iterations the solver may perform
     * 
     *  @return a solver instance whose initial guess is obtained by diagonalizing the core Hamiltonian
     */
    static DIISRHFSCFSolver<ExpansionScalar> WithCoreGuess(const RSpinorBasis<ExpansionScalar, GTOShell>& spinor_basis, const SQHamiltonian<ExpansionScalar>& sq_hamiltonian, const size_t N, const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6, const double threshold=1.0e-08, const size_t maximum_number_of_iterations = 128) {

        const auto C_initial = BaseRHFSCFSolver<ExpansionScalar>::calculateCoreGuess(sq_hamiltonian.core(), spinor_basis.overlap());

        return DIISRHFSCFSolver<ExpansionScalar>(C_initial, N, S, sq_hamiltonian, minimum_subspace_dimension, maximum_subspace_dimension, threshold, maximum_number_of_iterations);
    }


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return if the algorithm is considered to be converged
     */
    bool isConverged() override {
        if (this->numberOfIterations() == 0) {
            return false;
        } else {
            return ((this->D_current - this->D_previous).norm() < this->threshold);
        }
    }


    /**
     *  @return a new iterate to be used in the next iteration if the DIIS acceleration has not been switched on yet
     */
    // Iterate regularNewIterate() override {

    //     const auto F = QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(this->D_current, this->sq_hamiltonian);  // the Fock matrix constructed from the current coefficient matrix

    //     // Solve the generalized eigenvalue problem for the Fock matrix to get a new iteration of the coefficient matrix
    //     using MatrixType = Eigen::Matrix<ExpansionScalar, Eigen::Dynamic, Eigen::Dynamic>;
    //     Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (F.parameters(), this->S);
    //     const auto& C = generalized_eigensolver.eigenvectors();

    //     // Update the previous and current density matrices
    //     this->D_previous = this->D_current;
    //     this->D_current = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(C, this->N);

    //     return C;
    // }


    /**
     *  @return a new iterate according to the DIIS algorithm
     */
    Iterate updateIterate() override {

        // Start by creating a regular new iterate and calculating the associated error
        const auto F = QCModel::RHF<ExpansionScalar>::calculateScalarBasisFockMatrix(this->D_current, this->sq_hamiltonian);  // the Fock matrix constructed from the current coefficient matrix

        // Solve the generalized eigenvalue problem for the Fock matrix to get a new iteration of the coefficient matrix
        using MatrixType = Eigen::Matrix<ExpansionScalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixType> generalized_eigensolver (F.parameters(), this->S);
        const auto& C = generalized_eigensolver.eigenvectors();



        




        // Update the previous and current density matrices
        this->D_previous = this->D_current;
        this->D_current = QCModel::RHF<ExpansionScalar>::calculateScalarBasis1RDM(C, this->N);





        Error error = this->calculateError(iterate);
        this->error_deque.emplace(error);


        // Apply the DIIS acceleration to produce a better iterate if the current subspace dimension is large enough
        if (this->subspaceDimension() >= this->minimum_subspace_dimension) {

            const auto diis_coefficients = this->calculateDIISCoefficients();
            Iterate diis_iterate;  // the DIIS-accelerated iterate
            for (size_t i = 0; i < this->subspaceDimension(); i++) {
                diis_iterate += diis_coefficients(i) * this->iterate_deque[i];
            }
            iterate = diis_iterate;  // the accelerated iterate should be used in the next iteration, and it should also be added to the current subspace
        }
        this->iterate_deque.emplace(iterate);


        // Discard the oldest iterate and errors if the subspace becomes too large
        if (this->subspaceDimension() > this->maximum_subspace_dimension) {
            this->iterate_deque.pop_front();
            this->error_deque.pop_front();
        }
        return iterate;
    }




};


}  // namespace GQCP
